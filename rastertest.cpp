// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <fstream>
#include <cstdio>
#include <jpeglib.h>
#include "vecgen.hpp"
#include "raster.hpp"
#include "framebuf.hpp"
#include "texture.hpp"
#include "test.hpp"

const int tri1_w = 300, tri1_h = 300;

// Minimal libjpeg loader, jpg().load(buf,len,w,h) returns an array of
// w*h xRGB ints in natural (i.e. top-down) order.  Free with
// delete[].  Reuse the object for multiple files for a modest
// performance increase.
struct jpg {
    jpeg_source_mgr sm; // "superclass" must be first
    jpeg_decompress_struct jdec;
    jpeg_error_mgr errmgr;
    unsigned char* buf;
    int len;

    jpg() {
        memset(this, 0, sizeof(*this));
        sm.init_source = init_cb;
        sm.fill_input_buffer = fill_cb;
        sm.skip_input_data = skip_cb;
        jdec.err = jpeg_std_error(&errmgr);
        jpeg_create_decompress(&jdec);
        jdec.src = &sm;
    }

    ~jpg() { jpeg_destroy_decompress(&jdec); }

    unsigned int *load(unsigned char *mem, int bytes, int &w, int &h) {
        buf = mem;
        len = bytes;
        jpeg_read_header(&jdec, 1);
        jpeg_start_decompress(&jdec);
        w = jdec.output_width;
        h = jdec.output_height;
        if(jdec.output_components != 3)
            return 0;
        std::vector<unsigned char> line(3*w);
        unsigned int *img = new unsigned int[w*h];
        for(int y=0; y<h; y++) {
            unsigned char *sl = &line[0];
            jpeg_read_scanlines(&jdec, &sl, 1);
            for(int x=0; x<w; x++) {
                unsigned char *p = &line[x*3];
                img[w*y+x] = (p[0]<<16) | (p[1]<<8) | p[2];
            }
        }
        return img;
    }

    static void init_cb(j_decompress_ptr) {}

    static boolean fill_cb(j_decompress_ptr jdec) {
        jpg* me = (jpg*)jdec->src;
        jdec->src->next_input_byte = me->buf;
        jdec->src->bytes_in_buffer = me->len;
        return TRUE;
    }

    static void skip_cb(j_decompress_ptr jdec, long n) {
        jdec->src->next_input_byte += n;
        jdec->src->bytes_in_buffer -= n;
        if(jdec->src->bytes_in_buffer <= 0) {
            jdec->src->bytes_in_buffer = 0;
            jdec->src->next_input_byte = 0;
        }
    }
};

TEST(fb1) {
    framebuf fb(1, 1);

    unsigned int pixel;
    void *mems[] = { &pixel };
    const int memidx = 0;

    fb.set_color_mem(memidx);

    vecgen vg(optstring().c_str());

    vr color[] = { vg.input(0), vg.input(1), vg.input(2), vg.input(3) };
    float front_facing=0, x=0, y=0, z=0;
    fb.write(&vg, x, y, z, front_facing, color);

    float in[] = { 0.5, 0, 0.5, 1 };
    vecgen_run(vg, 0, in, 0, mems);

    return pixel == 0x007f007f;
}

static void dump_tri(const char *name, int w, int h, int *img)
{
    FILE* out = fopen(name, "wb");
    fprintf(out, "P6\n%d %d\n255\n", w, h);
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
            int rgb = img[y*w+x];
            fputc((rgb >> 16) & 0xff, out);
            fputc((rgb >> 8) & 0xff, out);
            fputc(rgb & 0xff, out);
        }
    }
    fclose(out);
}

struct tri1_t { vecgen *vg; framebuf *fb; };

void tri1_frag(vr sx, vr sy, vr *attdxdy, int natt, void *arg)
{
    (void)natt;
    tri1_t *t1 = (tri1_t*)arg;
    vr color[] = { attdxdy[3*0], attdxdy[3*1], attdxdy[3*2] };
    t1->fb->write(t1->vg, sx, sy, 0, 1, color);
}

TEST(tri1) {
    vr a[] = {  0,  1, 0.0, 1 };
    vr b[] = {  1, -1, 0.7, 1 };
    vr c[] = { -1,  0, 0.1, 1 };

    vr aatt[] = {1,1,0}, batt[] = {1,0,1}, catt[] = {0,1,1};

    // xRGB order, passed in as memidx zero
    unsigned int *tri1_img = new unsigned int[tri1_w*tri1_h];
    memset(tri1_img, 0, tri1_w*tri1_h*sizeof(int));
    void *mems[] = { tri1_img };

    vecgen vg(optstring().c_str());

    framebuf fb(tri1_w, tri1_h);
    fb.set_color_mem(0);
    tri1_t t1 = { &vg, &fb };

    rasterize_tri(&vg, 0, 0, tri1_w, tri1_h, a, b, c, 3, aatt, batt, catt,
                  tri1_frag, &t1);

    // Give the code a 16 slot output register set for debugging
    union { float f; int i; } outs[16];
    for(int i=0; i<16; i++)
        outs[i].i = 0x94beaaab;

    vecgen_run(vg, 0, 0, &outs[0].f, mems);

    // Log the values that got written
    for(int i=0; i<16; i++)
        if(outs[i].i != (int)0x94beaaab)
            printf("O%2.2d = 0x%8.8x (dec %d) (float %g)\n", i, outs[i].i, outs[i].i, outs[i].f);

    dump_tri("tri1.ppm", tri1_w, tri1_h, (int*)tri1_img);
    delete[] tri1_img;
    return true;
}

struct textri_t { vecgen *vg; framebuf *fb; texture *tex; };

void textri_frag(vr sx, vr sy, vr *attdxdy, int natt, void *arg)
{
    (void)natt;
    textri_t *tt = (textri_t*)arg;
    vr tc[]  = { attdxdy[0], attdxdy[3] };
    vr ddx[] = { attdxdy[1], attdxdy[4] };
    vr ddy[] = { attdxdy[2], attdxdy[5] };
    vr color[4];
    tt->tex->sample(tt->vg, 1, tc, ddx, ddy, color);
    swap(color[0], color[2]); // Framebuffer is 0xAARRGGBB, texture comes out {R,G,B,A}
    tt->fb->write(tt->vg, sx, sy, 0, 1, color);
}

TEST(textri) {
    std::ifstream f("textri.jpg");
    std::vector<char> v((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    int wh[2]; // extents array for texture
    int &w = wh[0], &h = wh[1]; // shorthand
    unsigned int *img = jpg().load((unsigned char*)&v[0], v.size(), w, h);

    texture tex(2, 4, texture::UBYTE);
    tex.set_mip_mode(texture::TRILINEAR);
    tex.set_dims(true, false, texture::EDGE);

    // Box filter (not quite right for NPOT: it produces a half-pixel
    // shift)
    for(int mip=0; w > 1 || h > 1; mip++) {
        int w2 = max(w/2, 1), h2 = max(h/2, 1);
        tex.set_mip(mip, wh, (char*)img);
        unsigned int *img2 = new unsigned int[w2 * h2];
        for(int y=0; y<h2; y++) {
            for(int x=0; x<w2; x++) {
                int p[] = { 2*w*y     + 2*x,  2*w*y     + (2*x+1),
                            2*w*(y+1) + 2*x,  2*w*(y+1) + (2*x+1) };
                int sums[4] = {};
                for(int c=0; c<4; c++)
                    for(int i=0; i<4; i++)
                        sums[c] += (img[p[i]] >> (c * 8)) & 0xff;
                for(int c=0; c<4; c++)
                    sums[c] = (sums[c] / 4) + !!(sums[c] & 2);
                img2[y*w2+x] = (sums[3]<<24)|(sums[2]<<16)|(sums[1]<<8)|sums[0];
            }
        }
        w = max(w2, 1); h = max(h2, 1);
        delete[] img;
        img = img2;
    }

    framebuf fb(tri1_w, tri1_h);
    unsigned int *fbimg = new unsigned int[tri1_w*tri1_h];
    memset(fbimg, 0, tri1_w*tri1_h*sizeof(int));
    fb.set_color_mem(0);

    void *mems[] = { fbimg, tex.get_mem() };

    vr a[] = {  0,  1, 0, 1 };
    vr b[] = {  1, -1, 0, 1 };
    vr c[] = { -1,  0, 0, 1 };
    vr aatt[] = {0.5,0}, batt[] = {1,1}, catt[] = {0,0.5};

    vecgen vg(optstring().c_str());
    textri_t tt = { &vg, &fb, &tex };
    rasterize_tri(&vg, 0, 0, tri1_w, tri1_h, a, b, c, 2, aatt, batt, catt,
                  textri_frag, &tt);

    vecgen_run(vg, 0, 0, 0, mems);

    dump_tri("textri.ppm", tri1_w, tri1_h, (int*)fbimg);

    return true;
}
