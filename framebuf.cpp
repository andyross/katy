// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include "vecgen.hpp"
#include "util.hpp"
#include "framebuf.hpp"

// FIXME: the stencil operations, with all their masking, are a mess
// when done with floating point math.  These will really benefit from
// integer ops in AVX2, but need to be factored out appropriately
// first...

const int ST_BITS = 8; // Only support 8 bit stencil right now

framebuf::framebuf(int height, int width)
{
    memset(this, 0, sizeof(*this));
    h = height;
    w = width;
    fb_memidx = db_memidx = st_memidx = -1;
}

void framebuf::set_depth(bool enabled, cmp_mode mode)
{
    depth_enabled = enabled;
    depth_mode = mode;
}

void framebuf::set_blend(bool color_channel, blend_eq eq, blend_fn src, blend_fn dst)
{
    if(color_channel) {
        blend_eq_c = eq;
        bf_srcc = src;
        bf_dstc = dst;
    } else {
        blend_eq_a = eq;
        bf_srca = src;
        bf_dsta = dst;
    }
}

void framebuf::set_constant_color(vr *const_color)
{
    for(int i=0; i<4; i++)
        blend_const[i] = const_color[i];
}

void framebuf::set_stencil(bool front, cmp_mode test, int ref, int mask,
                     stencil_op fail, stencil_op zfail, stencil_op zpass)
{
    stencil_mode &sm = front ? stencil_front : stencil_back;
    sm.mode = test;
    sm.ref = ref;
    sm.mask = mask;
    sm.fail = fail;
    sm.zfail = zfail;
    sm.zpass = zpass;
}

void framebuf::set_color_format(fb_format fmt, bool has_alpha)
{
    fb_fmt = fmt;
    fb_has_alpha = has_alpha;
}

void framebuf::set_depth_format(db_format fmt) { db_fmt = fmt; }
void framebuf::set_color_mem(int memidx) { fb_memidx = memidx; }
void framebuf::set_depth_mem(int memidx) { db_memidx = memidx; }
void framebuf::set_stencil_mem(int memidx) { st_memidx = memidx; }

bool framebuf::need_stencil()
{
    for(int i=0; i<2; i++) {
        stencil_mode *s = i ? &stencil_front : &stencil_back;
        if(s->mode != ALWAYS || s->fail != KEEP ||
           s->zfail != KEEP  || s->zpass != KEEP)
            return true;
    }
    return false;
}

vr framebuf::eval_cmp(cmp_mode m, vr a, vr b)
{
    switch(m) {
    case ALWAYS: return 1;
    case NEVER:  return 0;
    case LT:     return a < b;
    case LTE:    return a <= b;
    case GT:     return a > b;
    case GTE:    return a >= b;
    case EQ:     return a == b;
    case NEQ:    return a != b;
    };
    return 0;
}

void framebuf::load_depth_stencil(vr off)
{
    vr word = 0;
    memsz sz = db_fmt == INT16 ? U16 : U32;
    if(depth_enabled || (need_stencil() && db_fmt == PACKED8_24))
        word = load(sz, db_memidx, f2i(off));

    if(db_fmt == INT16) {
        db_tmp = i2f(word) * (1.0/0xffff);
    } else if(db_fmt == ZFLOAT) {
        db_tmp = word;
    } else {
        db_tmp = i2f(word & imm_i(0x00ffffff)) * (1.0/0xffffff);
        st_tmp = i2f(word & imm_i(0xff000000)) * (1.0/(1<<24));
    }

    if(need_stencil() && db_fmt != PACKED8_24)
        st_tmp = load(U8, st_memidx, f2i(off));

    ds_loaded = true;
}

vr framebuf::stencil_test(vr off, stencil_mode *st)
{
    if(!ds_loaded)
        load_depth_stencil(off);
    vr a, b;
    int fullmask = (1 << ST_BITS) - 1;
    if((st->mask & fullmask) == fullmask) {
        // All bits are in mask, we can skip it
        a = st->ref;
        b = st_tmp;
    } else {
        a = vg->imm(st->ref & st->mask);
        b = i2f(f2i(st_tmp) & imm_i(st->mask));
    }
    return eval_cmp(st->mode, a, b);
}

vr framebuf::get_depth(vr off)
{
    if(!ds_loaded)
        load_depth_stencil(off);
    return db_tmp;
}

vr framebuf::compute_blend(blend_fn &bf, int i, vr *dst, vr *src, vr *src1)
{
    if(bf.mode == MUL_ZERO)
        return 0;
    if(bf.mode == MUL_ONE)
        return src[i];

    vr f;
    i = bf.fn_use_color ? i : 3;
    switch(bf.fn_src) {
    case SRC:   f = src[i];         break;
    case SRC1:  f = src1[i];        break;
    case DST:   f = dst[i];         break;
    case CONST: f = blend_const[i]; break;
    }

    if(bf.fn_reverse)
        f = vr(1) - f;

    return f * src[i];
}

void framebuf::do_blend(vr off, vr *color, vr *color1)
{
    vr out[4];
    if(!blend_enabled) {
        for(int i=0; i<3; i++)
            out[i] = color[i];
        if(fb_has_alpha)
            out[3] = color[3];
    } else {
        // Read it
        vr fb[4];
        if(fb_fmt == ARGB8888) {
            vr word = load(U32, fb_memidx, f2i(off));
            fb[0] = (word & imm_i(0x00ff0000)) * (1.0/(255.*(1<<16)));
            fb[1] = (word & imm_i(0x0000ff00)) * (1.0/(255.*(1<<8)));
            fb[2] = (word & imm_i(0x000000ff)) * (1.0/255.);
            if(fb_has_alpha)
                fb[3] = (word & imm_i(0xff000000)) * (1.0/(255.*(1<<24)));
        } else if(fb_fmt == RGB565) {
            vr word = load(U16, fb_memidx, f2i(off));
            fb[0] = (word & imm_i(31<<11)) * (1.0/(31*(1<<11)));
            fb[1] = (word & imm_i(63<<5))  * (1.0/(63*(1<<5)));
            fb[2] = (word & imm_i(31))     * (1.0/31);
        } else if(fb_fmt == ARGB1555){
            vr word = load(U16, fb_memidx, f2i(off));
            fb[0] = (word & imm_i(31<<10)) * (1.0/(31*(1<<10)));
            fb[1] = (word & imm_i(31<<5))  * (1.0/(31*(1<<5)));
            fb[2] = (word & imm_i(31))     * (1.0/31);
            if(fb_has_alpha)
                fb[3] = mux(word & imm_i(1<<15), 1, 0);
        } else {
            vr base = off * (fb_has_alpha ? 4 : 3);
            for(int i=0; i<3; i++)
                fb[i] = load(U32, fb_memidx, f2i(base+i));
            if(fb_has_alpha)
                fb[3] = load(U32, fb_memidx, f2i(base+3));
        }

        if(!fb_has_alpha)
            fb[3] = 1;

        // Composite
        vr src[4], dst[4];
        for(int i=0; i<3; i++) {
            src[i] = compute_blend(bf_srcc, i, fb, color, color1);
            dst[i] = compute_blend(bf_dstc, i, fb, color, color1);
        }
        src[3] = compute_blend(bf_srca, 3, fb, color, color1);
        dst[3] = compute_blend(bf_dsta, 3, fb, color, color1);

        for(int i=0; i<4; i++) {
            switch(i != 3 ? blend_eq_c : blend_eq_a) {
            case SRC_PLUS_DST:  out[i] = src[i] + dst[i];     break;
            case SRC_MINUS_DST: out[i] = src[i] - dst[i];     break;
            case DST_MINUS_SRC: out[i] = dst[i] - src[i];     break;
            case MIN:           out[i] = min(src[i], dst[i]); break;
            case MAX:           out[i] = max(src[i], dst[i]); break;
            }
        }
    }

    // Write it out
    if(fb_fmt == ARGB8888) {
        vr r = vr(ulp_less(256) * (1<<16)) * max(min(out[0], 1), 0);
        vr g = vr(ulp_less(256) * (1<<8))  * max(min(out[1], 1), 0);
        vr b = vr(ulp_less(256))           * max(min(out[2], 1), 0);
        vr word = ((f2i(r) & imm_i(0x00ff0000)) |
                   (f2i(g) & imm_i(0x0000ff00)) |
                   (f2i(b) & imm_i(0x000000ff)));
        if(fb_has_alpha) {
            vr a = vr(ulp_less(256) * (1<<24)) * max(min(out[3], 1), 0);
            word |= f2i(a) & imm_i(0xff000000);
        }
        store(U32, fb_memidx, f2i(off), word);
    } else if(fb_fmt == RGB565) {
        vr r = vr(ulp_less(32) * (1<<11)) * max(min(out[0], 1), 0);
        vr g = vr(ulp_less(64) * (1<<5))  * max(min(out[1], 1), 0);
        vr b = vr(ulp_less(32))           * max(min(out[2], 1), 0);
        vr word = ((f2i(r) & imm_i(31<<11)) |
                   (f2i(g) & imm_i(63<<5)) |
                   (f2i(b) & imm_i(31)));
        store(U16, fb_memidx, f2i(off), word);
    } else if(fb_fmt == ARGB1555) {
        vr r = vr(ulp_less(32) * (1<<10)) * max(min(out[0], 1), 0);
        vr g = vr(ulp_less(32) * (1<<5))  * max(min(out[1], 1), 0);
        vr b = vr(ulp_less(32))           * max(min(out[2], 1), 0);
        vr word = ((f2i(r) & imm_i(31<<10)) |
                   (f2i(g) & imm_i(31<<5)) |
                   (f2i(b) & imm_i(31)));
        if(fb_has_alpha)
            word |= mux(out[3] >= 0.5, imm_i(1<<15), imm_i(0));
        store(U16, fb_memidx, f2i(off), word);
    } else if(fb_fmt == FLOAT) {
        vr base = off * (fb_has_alpha ? 4 : 3);
        for(int i=0; i<3; i++)
            store(U32, fb_memidx, f2i(base+i), out[i]);
        if(fb_has_alpha)
            store(U32, fb_memidx, f2i(base+3), out[3]);
    }
}

vr framebuf::eval_stencil(stencil_mode &sm, stencil_op op)
{
    vr mask = imm_i((1 << ST_BITS) - 1);
    vr s;
    switch(op) {
    case KEEP:      s = st_tmp; break;
    case SET_ZERO:  s = 0; break;
    case REPLACE:   s = sm.ref; break;
    case INCR:      s = f2i(min(i2f(st_tmp)+1, mask)); break;
    case INCR_WRAP: s = f2i(i2f(st_tmp)+1) & mask; break;
    case DECR:      s = f2i(max(i2f(st_tmp)-1, 0)); break;
    case DECR_WRAP: s = f2i(i2f(st_tmp)-1) & mask; break;
    case INVERT:    s = st_tmp ^ mask; break;
    }
    if(sm.mask != ((1<<ST_BITS)-1))
        s &= sm.mask;
    return s;
}

// Stencil update.  This is a big mess.  In principle, we need to
// compute six (!) different results (front/back cross the three
// result types) and mux between them at runtime.  In practice, the
// backface cases can often be eliminted at compile time, and many of
// the cases are going to be degenerate (e.g. KEEP)
vr framebuf::compute_stencil(const vr &front, vr sok, vr zok)
{
    vr st_fail[2], st_zfail[2], st_zpass[2];
    for(int i=0; i<2; i++) {
        stencil_mode &sm = i ? stencil_back : stencil_front;
        st_fail[i] = eval_stencil(sm, sm.fail);
        st_zfail[i] = eval_stencil(sm, sm.zfail);
        st_zpass[i] = eval_stencil(sm, sm.zpass);
    }
    if(!fixed_front) {
        st_fail[0] = mux(front, st_fail[0], st_fail[1]);
        st_zfail[0] = mux(front, st_zfail[0], st_zfail[1]);
        st_zpass[0] = mux(front, st_zpass[0], st_zpass[1]);
    }
    if(depth_enabled)
        return mux(sok, mux(zok, st_zpass[0], st_zfail[0]), st_fail[0]);
    else
        return mux(sok, st_zpass[0], st_fail[0]);
}

void framebuf::write(vecgen *vg_in, vr x, vr y, vr z, const vr &front, vr *color, vr *color1)
{
    vg = vg_in;
    new (&db_tmp) vr(vg->scratch()); // placement new hackery to get fresh scratch regs
    new (&st_tmp) vr(vg->scratch());
    ds_loaded = false;

    // Sanify inputs, might be immediates
    x.setvg(vg); y.setvg(vg); z.setvg(vg);
    for(int i=0; i<4; i++) {
        color[i].setvg(vg);
        if(color1)
            color1[i].setvg(vg);
    }

    // FIXME: this is sorta wrong.  It's possible to set an immediate
    // *false* value for front, meaning that culling is set to send
    // only backfacing polygons.  Which would be really weird, but
    // legal.  And it will break this logic.
    fixed_front = vg->is_immediate(front);

    vr off = y * w + x;

    // Stencil test
    vr sok = vg->imm(1);
    if(need_stencil()) {
        if(fixed_front) {
            sok = stencil_test(off, &stencil_front);
        } else {
            VRIF  (*vg, front) sok = stencil_test(off, &stencil_front);
            VRELSE(*vg)        sok = stencil_test(off, &stencil_back);
        }
    }

    // Depth test
    vr zok = vg->imm(1);
    if(depth_enabled)
        VRIF(*vg, sok)
            zok = eval_cmp(depth_mode, z, get_depth(off));

    // Compute new stencil value
    vr st = 0;
    if(need_stencil())
        st = compute_stencil(front, sok, zok);

    // Z/stencil write
    if(depth_enabled || need_stencil()) {
        vr zword = z;
        if(db_fmt == PACKED8_24) {
            zword = max(0, min(z, 1)) * 0xffffff;
            if(need_stencil())
                zword |= f2i(i2f(st) * (1<<24));
        } else if(db_fmt == INT16) {
            zword = max(0, min(z, 1)) * 0xffff;
            store(U16, db_memidx, f2i(off), zword);
        }
        store(db_fmt == INT16 ? U16 : U32, db_memidx, f2i(off), zword);
        if(need_stencil() && db_fmt != PACKED8_24)
            store(U8, st_memidx, f2i(off), st);
    }

    // Framebuffer blend & write
    VRIF(*vg, sok && zok)
        do_blend(off, color, color1);

    vg = 0;
}
