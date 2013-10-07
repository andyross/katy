// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _TEXTURE_H
#define _TEXTURE_H

#include "vecgen.hpp"

class texture {
public:
    enum component_type { UBYTE, USHORT, INT, FLOAT };
    enum mip_mode { NOMIP, NEAREST, TRILINEAR };
    enum repeat_mode { REPEAT, EDGE, BORDER };

    // Usage note: components must always pack into a integer number
    // of the nearest LOADable type that fits.  So for example
    // 2-component UBYTE textures are fine but "RGB" is disallowed,
    // because you need 4 to fill the load.  Use a 4-component texture
    // and just ignore the alpha channel, the optimizer will take it
    // out.
    texture(unsigned int dim, unsigned int ncomponents, component_type c);
    ~texture();
    void set_mip_mode(mip_mode m) { mip = m; }
    void set_dim(unsigned int d, bool linear, bool mirror, repeat_mode wrap);
    void set_dims(bool linera, bool mirror, repeat_mode wrap);
    void set_mip(unsigned int level, int *extents, char* mem);
    void set_max_aniso_samples(int n) { max_aniso_samples = n; }
    void set_lod_clamp(int min, int max) { min_mip = min; max_mip = max; }
    void set_border_color(float *c);
    void set_lod_premul(float p) { mip_premul = p; }

    // The texture stores itself in a single memory block (extents
    // metadata first, followed by the texel data).  The user needs to
    // fetch this pointer (to be placed in the mem[] array of the
    // generated function), and pass its index to load().  Note that
    // the pointer will change at each call to set_mip().
    void *get_mem() { return vg_mem; }

    bool validate();
    void sample(vecgen *vg, int memidx, vr *coords, vr *ddx, vr *ddy, vr *out);

private:
    struct dimension {
        bool linear;
        bool mirror;
        repeat_mode wrap;
    };

    struct miprec {
        std::vector<int> extents; // {x/width, y/height, z/depth, ... }
        unsigned int start;           // index into memory array (load units)
    };

    // Input parameters:
    unsigned int dim;
    unsigned int n_components;
    component_type ctype;
    mip_mode mip;
    std::vector<dimension> dims; // indexed by dimension
    std::vector<miprec> mips; // indexed by mip level

    // FIXME: these are good candidates for specification at runtime
    // as vr's, not just as immediates at codegen time.
    int max_aniso_samples;
    int min_mip, max_mip;       // LOD clamping
    std::vector<float> border_color; // indexed by n_components
    float mip_premul;           // LOD bias

    char *vg_mem;  // raw memory array
    vecgen *vg;    // cached load() argument
    int vg_memidx; // cached load() argument

    void load_mipped(vr *coords, vr *ddx, vr *ddy, vr *out);
    void load_2mips(vr mip0, vr mip1, vr frac, vr *coords, vr *out);
    void load_1mip(vr &mip, vr *coords, vr *out);
    void load_simple(vr *coords, vr *offsets, vr start_idx, vr *extents, vr *out);
    vr read_mip(vr &mip, vr *extents);
    vr smart_recip(double max, vr a) { return max > 1024 ? vg->recip2(a) : recip(a); }
    vr abs(vr a) { return mux(a < imm(0), imm(-1)*a, a); }
    static int compsz(component_type c);
};

#endif // _TEXTURE_H
