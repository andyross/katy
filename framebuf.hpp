// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _FRAMEBUF_HPP
#define _FRAMEBUF_HPP

// Hopefully OpenGL 4.2-compilant framebuffer code generation
class framebuf {
public:
    enum cmp_mode { ALWAYS, NEVER, LT, LTE, GT, GTE, EQ, NEQ };
    enum stencil_op { KEEP, SET_ZERO, REPLACE, INCR, INCR_WRAP, DECR, DECR_WRAP, INVERT };
    enum blend_eq { SRC_PLUS_DST, SRC_MINUS_DST, DST_MINUS_SRC, MIN, MAX };
    enum blend_fn_mode { MUL_ZERO, MUL_ONE, FUNC };
    enum blend_fn_src { SRC, DST, CONST, SRC1 };
    enum fb_format { ARGB8888, RGB565, ARGB1555, FLOAT };
    enum db_format { PACKED8_24, INT16, ZFLOAT };

    struct blend_fn {
        blend_fn_mode mode;
        blend_fn_src fn_src;
        bool fn_use_color; // false == alpha
        bool fn_reverse; // true == ONE_MINUS_x
    };

    framebuf(int h, int w);

    void set_depth(bool enabled, cmp_mode mode);
    void set_blend(bool color_channel, blend_eq eq, blend_fn src, blend_fn dst);
    void set_constant_color(vr *const_color);

    // Set stencil parameters (KEEP/ALWAYS for all values means "disabled")
    // FIXME: should ref/mask be vr's, to allow for runtime setting?
    void set_stencil(bool front, cmp_mode test, int ref, int mask,
                     stencil_op fail, stencil_op zfail, stencil_op zpass);

    void set_color_format(fb_format fmt, bool has_alpha);
    void set_depth_format(db_format fmt);
    void set_color_mem(int memidx);
    void set_depth_mem(int memidx);
    void set_stencil_mem(int memidx);

    void write(vecgen *vg_in, vr x, vr y, vr z, const vr &front, vr *color, vr *color1=0);

private:
    struct stencil_mode {
        cmp_mode mode;
        int ref;
        int mask;
        stencil_op fail, zfail, zpass;
    };

    void load_depth_stencil(vr off);
    vr stencil_test(vr off, stencil_mode *st);
    vr eval_cmp(cmp_mode m, vr a, vr b);
    vr get_depth(vr off);
    bool need_stencil();
    vr compute_stencil(const vr &front, vr sok, vr zok);
    vr compute_blend(blend_fn &bf, int i, vr *dst, vr *src, vr *src1);
    void do_blend(vr off, vr *color, vr *color1);
    vr eval_stencil(stencil_mode &sm, stencil_op op);

    int w, h;

    stencil_mode stencil_front, stencil_back;

    bool depth_enabled;
    cmp_mode depth_mode;

    bool blend_enabled;
    blend_eq blend_eq_a, blend_eq_c;
    blend_fn bf_srca, bf_srcc, bf_dsta, bf_dstc;
    vr blend_const[4];

    fb_format fb_fmt;
    bool fb_has_alpha;
    db_format db_fmt;
    int fb_memidx, db_memidx, st_memidx;

    vecgen *vg;
    vr db_tmp, st_tmp;
    bool ds_loaded;
    bool fixed_front; // only front-facing, known at compile time
};

#endif // _FRAMEBUF_HPP
