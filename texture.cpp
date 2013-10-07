// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <cmath>
#include <vector>
#include <algorithm>
#include "util.hpp"
#include "texture.hpp"

using namespace std;

// FIXME: cube maps (and in particular wrap-to-adjacent-face)

// FIXME: work up the "CUSTOM" handling for things like texture
// compression (needs to take a single index and do N loads, returning
// an array of output colors, should have some checking
// vs. n_components. etc...)

// FIXME: OpenGL requires separate minification/magnification
// linearity, but this code is simply "linear or not". (note that i915
// apparently does not support NEAREST magnification at all, and it
// gets away with it...)

// FIXME: RGB332, RGB565, ARGB1555, ARGB4444, ARGB2:10:10:10, HALF
// (when supported), YUYV, "CUSTOM" with plugin functionality

// FIXME: support for putting the data in externally-managed memory
// (i.e. shared mappings for buffer objects)

texture::texture(unsigned int dim, unsigned int n, component_type c)
    : dim(dim), n_components(n), ctype(c), mip(NOMIP),
      max_aniso_samples(1), min_mip(-1), max_mip(-1), mip_premul(1),
      vg_mem(0)
{
    dimension d = { false, false, REPEAT };
    dims.assign(dim, d);
}

texture::~texture()
{
    delete[] vg_mem;
}

void texture::set_dim(unsigned int d, bool linear, bool mirror, repeat_mode wrap)
{
    dimension dim = { linear, mirror, wrap };
    dims[d] = dim;
}

void texture::set_dims(bool l, bool m, repeat_mode w)
{
    for(unsigned int i=0; i<dim; i++)
        set_dim(i, l, m, w);
}

int texture::compsz(component_type c)
{
    if(c == UBYTE) return 1;
    if(c == USHORT) return 2;
    return 4;
}

// This repacks the array every call.  Could be done much faster by
// presizing and appending, but texture initialiation isn't normally a
// performance path.
void texture::set_mip(unsigned int level, int *extents, char* mem)
{
    if(level >= mips.size())
        mips.resize(level+1);

    mips[level].extents.resize(dim);
    for(unsigned int i=0; i<dim; i++)
        mips[level].extents[i] = extents[i];

    int tsz = n_components * compsz(ctype); // bytes per texel
    int lu = tsz < 4 ? tsz : 4; // bytes per "load unit"
    int lsz = tsz / lu; // load units per texel

    // Pack the metadata
    vector<uif> meta;
    for(unsigned int i=0; i<mips.size(); i++) {
        uif u;
        if(dim == 2) {
            u.i = (mips[i].extents[1] << 16) | mips[i].extents[0];
            meta.push_back(u);
        } else {
            for(unsigned int d=0; d<dim; d++) {
                u.f = mips[i].extents[d];
                meta.push_back(u);
            }
        }
        meta.push_back(u); // make room for start index
    }

    // Size the array.  All sizes below in load units
    int metasz = meta.size() * sizeof(meta[0]) / lu;
    int sz = metasz;
    vector<int> starts, sizes;
    for(unsigned int i=0; i<mips.size(); i++) {
        int mipsz = 0;
        if(mips[i].extents.size() == dim) {
            mipsz = lsz;
            for(unsigned int d=0; d<dim; d++)
                mipsz *= mips[i].extents[d];
        }
        sizes.push_back(mipsz);
        starts.push_back(sz);
        meta[((dim < 3 ? 1 : dim)+1)*(i+1)-1].f = sz;
        sz += mipsz;
    }

    // Allocate and copy
    char *newmem = new char[(metasz+sz)*lu];
    memcpy(newmem, &meta[0], metasz*lu);
    for(unsigned int i=0; i<mips.size(); i++) {
        if(mips[i].extents.size() == dim) {
            char *src = i == level ? mem : &vg_mem[lu*mips[i].start];
            memcpy(newmem + lu*starts[i], src, lu*sizes[i]);
        }
        mips[i].start = starts[i];
    }
    delete[] vg_mem;
    vg_mem = newmem;
}

void texture::set_border_color(float *c)
{
    border_color.clear();
    for(unsigned int i=0; i<dim; i++)
        border_color.push_back(c[i]);
}

bool texture::validate()
{
    if(max_aniso_samples > 1 && dim != 2)
        return false;

    // Check the mipmap array that it shrinks by half each step.
    // These use the ARB NPOT semantics: each mip level dimension must
    // be the floor of half its parent, clamped to a minimum of 1.
    if(mips.size() < 1) return false;
    for(unsigned int i=1; i<mips.size(); i++) {
        if(mips[i].extents.size() != dim)
            return false;
        for(unsigned int d=0; d<dim; d++)
            if(mips[i].extents[d] != max(mips[i-1].extents[d]/2, 1))
                return false;
    }

    // We use floating point math to compute indexes, so that sets a
    // bound on how big the array can be.
    int sztotal = 0;
    for(unsigned int i=1; i<mips.size(); i++) {
        int sz = (n_components * compsz(ctype) + 3)/4;
        for(unsigned int d=0; d<dim; d++)
            sz *= mips[i].extents[d];
        sztotal += sz;
    }
    if(sztotal > 0xffffff)
        return false;

    if(min_mip < 0 || min_mip >= (int)mips.size()) min_mip = 0;
    if(max_mip < 0 || max_mip > min_mip) max_mip = mips.size()-1;
    return true;
}

static float pow2(int p)
{
    union { int i; float f; } u;
    u.i = (p+127) << 23;
    return u.f;
}

// At runtime, the mip metadata lives at the start of the texture
// buffer as an array of 32 bit floats.  There is a record of "dim+1"
// floats for each mip, in order.  For each mip: first comes an array
// of "dim" extent (width/height/etc...) values.  Then comes a start
// index indicating the location (in LOAD units) within the buffer
// where the texel data for the mip begins.
//
// EXCEPTION: for the special case of 2D textures, the extents are not
// separate floats but instead packed (little endian) as two 16 bit
// *integer* values, so that they can be loaded simultaneously.  3D+
// extents cannot be packed that way in the general case; they may be
// too large to fit. FIXME: detect small textures and bit-pack 3D+
// extents when they fit.
vr texture::read_mip(vr &mip, vr *extents)
{
    if(vg->is_immediate(mip)) {
        // Non-mipped loads technically use this API, but pass imm(0)
        // here, so handle specially for efficiency (eliding the two
        // LOAD ops is a big win).
        for(unsigned int i=0; i<dim; i++)
            extents[i] = imm(mips[0].extents[i]);
        return imm(mips[0].start);
    }

    int recsz = 1 + (dim > 2 ? dim : 1);
    vr recidx = mip * vg->imm(recsz);
    if(dim == 2) {
        vr tmp = load(U32, vg_memidx, f2i(recidx));
        extents[0] = i2f(tmp & imm_i(0xffff));
        extents[1] = i2f(tmp & imm_i(0xffff0000)) * pow2(-16);
    } else {
        for(unsigned int i=0; i<dim; i++)
            extents[i] = load(U32, vg_memidx, f2i(i ? recidx : recidx + i));
    }
    return load(U32, vg_memidx, f2i(recidx + (recsz - 1)));
}

// Generalized unfiltered texture load handles repeating, clamping to
// edge or border (and mirroring of each of those), and arbitrary
// output component counts and sizes.
void texture::load_simple(vr *coords, vr *offsets, vr start_idx, vr *extents, vr *out)
{
    // Fix up coordinates
    vector<vr> c(dim);
    c.assign(coords, coords+dim);
    bool use_border = false;
    vr border;
    border = vg->imm(0);
    for(unsigned int i=0; i<dim; i++) {
        // Annoying recip (and slow for big textures!).  Seems like
        // that should be fixable but wrap/clamp needs to be in
        // normalized space...
        c[i] += offsets[i] * smart_recip(mips[0].extents[i], extents[i]);
        use_border = use_border || dims[i].wrap == BORDER;
        if(dims[i].mirror)
            c[i] *= 0.5;
        if(dims[i].wrap == REPEAT)
            c[i] -= floor(c[i]);
        if(dims[i].wrap == BORDER)
            border = border || c[i] > ulp_less(1.0) || c[i] < 0;
        if(dims[i].mirror)
            c[i] = imm(1.0) - abs(c[i] - 2.0);
        c[i] = min(ulp_less(1.0), max(0, c[i]));
        c[i] = floor(c[i] * extents[i]);
    }

    // Compute the index.  By convention the first coordinate (x/u)
    // has stride 1, so walk backwards down dim.
    vr idx = c[dim-1];
    for(int i=dim-2; i>=0; i--)
        idx = idx * extents[i] + c[i];

    // Do the load...
    int bytes = compsz(ctype);
    int nl = max(n_components * bytes / 4, 1u);
    int lu = bytes * n_components / nl;
    memsz msz = lu==1 ? U8 : (lu==2 ? U16 : U32);
    idx = idx * nl + start_idx;
    vector<vr> loaded(nl);
    for(int i=0; i<nl; i++)
        loaded[i] = load(msz, vg_memidx, f2i(idx + i));

    if(bytes == 4) {
        for(unsigned int i=0; i<n_components; i++)
            out[i] = loaded[i];
    } else {
        // Goofball generic size vs. count logic.  But it works.  Note
        // that the subcomponent order is little endian.
        // FIXME: for byte loads, it ends up doing (load() & 0xff)
        // which is redundant.  Either fix here or add a case in the
        // optimier to recognize that mask for 8bit loads.
        int oi=0, mask=0, li=-1, shift=0;
        while(oi < (int)n_components) {
            if(!mask) {
                mask = (1<<(8*bytes))-1;
                shift = 0;
                li++;
            }
            out[oi++] = i2f(loaded[li] & imm_i(mask)) * pow2(-shift);
            mask <<= 8*bytes;
            shift += 8*bytes;
        }
    }

    // Emit conversion logic where needed (note only INT gets the
    // i2f(), because bytes and shorts do it above)
    for(unsigned int i=0; i<n_components; i++) {
        switch(ctype) {
        case UBYTE:  out[i] *= 1.0/0xff;                          break;
        case USHORT: out[i] *= 1.0/0xffff;                        break;
        case INT:    out[i] = i2f(out[i]) * ulp_less(pow2(-31));  break;
        case FLOAT: break;
        }
    }

    // ... or not.  (This always does the load for border pixels and
    // throws away the result, on the assumption that border pixels
    // are rare and that the extra work to mask out the load would be
    // a net harm.)
    if(use_border)
        for(unsigned int i=0; i<n_components; i++)
            out[i] = mux(border, border_color[i], out[i]);
}

static vr lerp(vr &a, vr &b, vr &frac)
{
    return fma(a, frac, b-a);
}

// Generalized linear filtering implementation works for all component
// counts and dimentionalities, and falls back to uninterpolated
// output (per-axis) if linear is not enabled.
void texture::load_1mip(vr &mip, vr *coords, vr *out)
{
    // FIXME: the read_mip() call gets needlessly repeated across
    // multiple load_2mips() calls when aniso is enabled.  Do
    // the call there and just pass down index and extents...
    vector<vr> extents(dim);
    vr idx = read_mip(mip, &extents[0]);

    vector<vr> coeffs;
    vector<vector<vr> > corners;
    corners.push_back(vector<vr>(dim)); // assign scratch registers
    float offset = dims[0].linear ? -0.5 : 0;
    corners.back().assign(dim, vg->imm(offset)); // assignment op, not constructor
    for(unsigned int d=0; d<dim; d++) {
        if(dims[d].linear) {
            // Duplicate all the corners, setting the Dth coordinate
            // to 1.0 in the new subspace.
            unsigned int n = corners.size();
            for(unsigned int i=0; i<n; i++) {
                corners.push_back(vector<vr>(dim));
                corners.back().assign(corners[i].begin(), corners[i].end());
                corners.back()[d] = imm(0.5);
            }
            vr diff = extents[d]*coords[d] + 0.5;
            diff -= floor(diff);
            coeffs.push_back(diff);
        }
    }

    // Do all the texture loads
    vector<vector<vr> > colors(corners.size());
    for(unsigned int i=0; i<corners.size(); i++) {
        colors[i].resize(n_components);
        load_simple(coords, &corners[i][0], idx, &extents[0], &colors[i][0]);
    }

    // Now "reduce" by interpolating along each axis, halving the size array each step
    for(/**/; coeffs.size(); coeffs.pop_back()) {
        vr &c = coeffs.back();
        int n = (coeffs.size()+1)/2;
        for(int i=0; i<n; i++)
            for(unsigned int j=0; j<n_components; j++)
                colors[i][j] = lerp(colors[i][j], colors[i+n][j], c);
    }
    for(unsigned int i=0; i<n_components; i++)
        out[i] = colors[0][i];
}

void texture::load_2mips(vr mip0, vr mip1, vr frac, vr *coords, vr *out)
{
    // Start with the first mip level...
    load_1mip(mip0, coords, out);
    if(mip != TRILINEAR)
        return;

    // ...and blend in the second.  Use the VRIF here to avoid doing
    // needless LOADs in the case where we have only one usable level
    // (mag filtering, lod clamping)
    VRIF(*vg, mip0 != mip1) {
        vector<vr> out1(n_components);
        load_1mip(mip1, coords, &out1[0]);
        for(unsigned int i=0; i<dim; i++)
            out[i] = lerp(out[i], out1[i], frac);
    }
}

void texture::load_mipped(vr *coords, vr *ddx, vr *ddy, vr *out)
{
    bool do_aniso = max_aniso_samples > 1 && dim == 2;

    // Texture axis derivatives. Compute du/dv/etc... as normalized
    // vectors in texture space (i.e. at magnitude 1.0 they're equal
    // in size to a texel).  Get the "scale" as the minimum of the
    // squared lengths of these vectors.
    vector<vr> dtdx(dim), dtdy(dim); // FIXME: no need for vectors here, these are temporaries
    vr scale, axis_lsq, axisx, axisy;
    vr ex = mips[0].extents[0], ey = mips[0].extents[1];
    for(unsigned int i=0; i<dim; i++) {
        dtdx[i] = ex * ddx[i];
        dtdy[i] = ey * ddy[i];
        vr lsq = dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i];
        if(do_aniso) {
            // compute the long-axis scale and save off the axis directions
            axis_lsq = i ? max(lsq, axis_lsq) : lsq;
            axisx = i ? mux(axis_lsq > lsq, axisx, ddx[i]) : ddx[i];
            axisy = i ? mux(axis_lsq > lsq, axisy, ddy[i]) : ddy[i];
        }
        scale = i ? min(scale, lsq) : lsq;
    }
    vr scale_sq = scale;

    // Precision muckery to compute sqrt(scale).  For mip selection,
    // the low precision variants are fine, as we only want the
    // logarithm.  The mantissa is used as the trilinear interpolation
    // constant, where 11 bits is OK for typical UBYTE color
    // components, but almost certainly not if the user is expecting
    // 16+ bits of precision.
    //
    // FIXME: no need for the recip, just deal in negative log space
    // (reverse the subtraction for mip0, frac becomes 1-frac, I
    // think...).
    if(ctype == UBYTE || mip != TRILINEAR)
        scale = recip(rsqrt(scale));
    else
        scale = vg->recip2(vg->rsqrt2(scale));

    // Apply bias as a multipler in texel space.  OpenGL wants it to
    // be an additive value in mip space, but we never see a true
    // continuous value there, using a trick to get the logarithm.
    if(mip_premul != 1)
        scale *= mip_premul;

    // When trilinear is enabled, we want an precise integer scale
    // (i.e. where the pixel size equals the texel size) to be the
    // point where mip choices change (from N-1,N to N,N+1).  When mip
    // choices are discrete, we want this to be the *center* of the
    // choice range.  Premultiply scale (by sqrt2, not 1.5, because
    // the interpolation is in log space!).
    if(mip == NEAREST)
        scale *= sqrt(2);

    vr mip0;
    mip0 = vg->imm(0);
    if(mip != NOMIP) {
        mip0 = i2f(scale & imm_i(0x7f800000)) * pow2(-23) - 127;
        mip0 = max(mip0, min_mip);
        mip0 = min(mip0, max_mip);
    }

    if(mip == NEAREST && !do_aniso) {
        load_1mip(mip0, coords, out);
        return;
    }

    // Interpolation coefficient.  Really we should interpolate the
    // logarithm of the mantissa, but log2(x) ~= (x-1) over the range
    // [1,2], so just use the fraction field from the IEEE number.
    vr frac = i2f(scale & imm_i(0x007fffff)) * pow2(-23);
    frac = mux(mip0 < min_mip, 0, frac);
    frac = mux(mip0 > max_mip, 1, frac);

    vr mip1 = mip0 + 1;
    mip1 = mux(mip1 > max_mip, max_mip, mip1);

    // Simple trilinear filtering can stop here
    if(!do_aniso) {
        load_2mips(mip0, mip1, frac, coords, out);
        return;
    }

    // Anisotropic filtering uses an iterated loop at runtime (because
    // most loads don't need all samples)
    vr samples = rsqrt(scale_sq*recip(axis_lsq)); // =~ sqrt(axis_lsq)/sqrt(scale)
    samples = min(max_aniso_samples, max(samples + 0.5, 1));

    // Interpolate aniso steps across +/-0.5 * axis
    vr isamp = recip(samples);
    vr stepx = axisx * isamp;
    vr stepy = axisy * isamp;
    vr c[2];
    c[0] = coords[0] + (stepx - axisx)*0.5;
    c[1] = coords[1] + (stepy - axisy)*0.5;

    for(unsigned int i=0; i<n_components; i++)
        out[i] = 0;
    vector<vr> tmp(n_components);
    VRLOOP(*vg) {
        load_2mips(mip0, mip1, frac, c, &tmp[0]);
        for(unsigned int i=0; i<n_components; i++)
            out[i] += tmp[i];
        c[0] += stepx;
        c[1] += stepy;
        samples -= 1;
        VRIF(*vg, samples <= 0)
            vg->break_loop();
    }
    for(unsigned int i=0; i<n_components; i++)
        out[i] *= isamp;
}

void texture::sample(vecgen *vgen, int mi, vr *coords, vr *ddx, vr *ddy, vr *out)
{
    if(validate()) {
        vg = vgen;
        vg_memidx = mi;
        for(unsigned int i=0; i<dim; i++) {
            coords[i].setvg(vg);
            if(ddx) ddx[i].setvg(vg);
            if(ddy) ddy[i].setvg(vg);
        }
        for(unsigned int i=0; i<n_components; i++)
            out[i].setvg(vg);
        if(mip == NOMIP && max_aniso_samples == 1) {
            vr imm0 = imm(0);
            load_1mip(imm0, coords, out);
        } else {
            load_mipped(coords, ddx, ddy, out);
        }
        vg = 0;
    }
}
