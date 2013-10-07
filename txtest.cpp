// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <cstdio>
#include <cmath>
#include <string>
#include "vecgen.hpp"
#include "texture.hpp"
#include "test.hpp"

// FIXME: things untested:
// + validate() failure cases
// + read_mip() for dim != 2 (i.e. the simple version)
// + MIRROR wrapping
// + mip_premul (i.e. lod bias) != 1.0

using namespace std;

// Rough equality test with special casing for zero
bool approx(int bits, double a, double b) {
    const double prec = 1./(1<<bits);
    return (a == 0 && b == 0) || (b != 0 && fabs(1 - a/b) < prec);
}

TEST(aniso)
{
    texture tex(2, 1, texture::INT);
    // == 1/8th, 7/8th
    unsigned int buf[] = { 0x10000000, 0x70000000, 0x10000000, 0x70000000 };
    int wh[] = { 2, 2 };
    tex.set_mip(0, wh, (char*)&buf[0]);
    tex.set_max_aniso_samples(16);

    void* mems[] = { tex.get_mem() };
    const int memidx = 0;

    vecgen vg(optstring().c_str());
    vr out0 = vg.output(0);
    vr coords[2] = { vg.input(0), vg.imm(0.25) };
    vr ddx[] = { vg.input(1),     vg.imm(0) };
    vr ddy[] = { vg.imm(0), vg.input(2) };
    tex.sample(&vg, memidx, coords, ddx, ddy, &out0);
    vg.codegen();

    float out, in[3];
    in[1] = 0.5;    // ddx
    in[2] = 0.5/16; // ddy
    const int samples = 23;
    float p0 = buf[0] / (float)0x80000000u;
    float p1 = buf[1] / (float)0x80000000u;
    bool ok = true;
    for(int i=0; i<samples; i++) {
        float f = i / (samples - 1.0);
        // Interpolate across one pixel in X
        in[0] = 0.25 + f/wh[0];
        vecgen_run_no_gen(vg, 0, in, &out, mems);

        float expect = p0 * (1-f) + p1 * f;
        ok = ok && approx(3, expect, out); // 3 bits for aniso_samples==16
    }
    return ok;
}

TEST(mipsel)
{
    texture tex(2, 2, texture::FLOAT);
    float mip0[] = { 123,0, 123,0, 123,0, 123,0 }; // mip 0 is color [123,0]
    float mip1[] = { 0,198 };                      // mip 1 is color [0,198]
    int wh0[] = { 2, 2 };
    int wh1[] = { 1, 1 };
    tex.set_mip(0, wh0, (char*)&mip0);
    tex.set_mip(1, wh1, (char*)&mip1);

    void* mems[] = { tex.get_mem() };
    const int memidx = 0;

    float ins[4], out[2];
    bool ok = true;

    // Subtle: trilinear interpolation output is precision-limited by
    // RSQ/RCP.  The code will do newton/raphson when using FLOAT
    // components, but it still doesn't come out perfect.
    const int bits = 22;

    for(int i=0; i<2; i++) {
        const bool linear = i == 1;
        tex.set_mip_mode(linear ? texture::TRILINEAR : texture::NEAREST);

        vecgen vg(optstring().c_str());
        vr texels[] = { vg.output(0), vg.output(1) };
        vr ddx[] = { vg.input(0), vg.input(1) }, ddy[] = { vg.input(2), vg.input(3) };
        vr coords[] = { vg.imm(0.5), vg.imm(0.5) };
        tex.sample(&vg, memidx, coords, ddx, ddy, texels);
        vg.codegen();

        // One-texel (0.5x0.5) size should pick mip level 0
        ins[0] = 0.5; ins[1] = 0;   // ddx
        ins[2] = 0;   ins[3] = 0.5; // ddy
        vecgen_run_no_gen(vg, 0, ins, out, mems);
        ok = ok && approx(bits, out[0], 123) && approx(bits, 123-out[1], 123);

        // Two-texture size should pick mip level 1
        ins[0] = 1; ins[1] = 0;   // ddx
        ins[2] = 0; ins[3] = 1; // ddy
        vecgen_run_no_gen(vg, 0, ins, out, mems);
        ok = ok && approx(bits, 198-out[0], 198) && approx(bits, out[1], 198);

        // A tiny bit less than midway, check that NEAREST uses mip0
        // and that TRILINEAR correctly interpolates.
        const int bits2 = 2;
        float mid = 0.5 * sqrt(2) * .9995;
        ins[0] = mid; ins[1] = 0;
        ins[2] = 0;   ins[3] = mid;
        vecgen_run_no_gen(vg, 0, ins, out, mems);
        if(linear) ok = ok && approx(bits2, out[0], 123/2.) && approx(bits2, out[1], 198/2.);
        else       ok = ok && out[0] == 123 && out[1] == 0;

        // A tiny bit more than midway, check that NEAREST uses mip1
        mid = 0.5 * sqrt(2) * 1.0005;
        ins[0] = mid; ins[1] = 0;
        ins[2] = 0;   ins[3] = mid;
        vecgen_run_no_gen(vg, 0, ins, out, mems);
        if(linear) ok = ok && approx(bits2, out[0], 123/2.) && approx(bits2, out[1], 198/2.);
        else       ok = ok && out[0] == 0 && out[1] == 198;
    }
    return ok;
}

TEST(minimipsel)
{
    texture tex(2, 2, texture::FLOAT);
    float mip0[] = { 123,0, 123,0, 123,0, 123,0 }; // mip 0 is color [123,0]
    float mip1[] = { 0,198 };                      // mip 1 is color [0,198]
    int wh0[] = { 2, 2 };
    int wh1[] = { 1, 1 };
    tex.set_mip(0, wh0, (char*)&mip0);
    tex.set_mip(1, wh1, (char*)&mip1);

    void* mems[] = { tex.get_mem() };
    const int memidx = 0;

    float ins[4], out[2];

    tex.set_mip_mode(texture::TRILINEAR);

    vecgen vg(optstring().c_str());
    vr texels[] = { vg.output(0), vg.output(1) };
    vr ddx[] = { vg.input(0), vg.input(1) }, ddy[] = { vg.input(2), vg.input(3) };
    vr coords[] = { vg.imm(0.5), vg.imm(0.5) };
    tex.sample(&vg, memidx, coords, ddx, ddy, texels);
    vg.codegen();

    // One-texel (0.5x0.5) size should pick mip level 0
    ins[0] = 0.5; ins[1] = 0;   // ddx
    ins[2] = 0;   ins[3] = 0.5; // ddy
    vecgen_run_no_gen(vg, 0, ins, out, mems);
    return approx(11, out[0], 123) && approx(11, 123-out[1], 123);
}

TEST(linear1d)
{
    const int samples = 5;

    texture tex(1, 1, texture::UBYTE);
    tex.set_dims(true, false, texture::REPEAT);
    unsigned char data[] = { 0, 0xff };
    int wh[] = { 2 };
    tex.set_mip(0, wh, (char*)&data);

    void* mems[] = { tex.get_mem() };
    const int memidx = 0;

    vecgen vg(optstring().c_str());
    vr texels[samples], coords[2];
    float expect[samples];
    for(int i=0; i<samples; i++) {
        float u = 0.25 + 0.5 * i/(samples-1.0);
        expect[i] = i/(samples-1.0);
        coords[0] = vg.imm(u);
        tex.sample(&vg, memidx, coords, 0, 0, &texels[i]);
        vg.output(i) = texels[i];
    }

    float out[samples];
    vecgen_run(vg, 0, 0, out, mems);
    bool ok = true;
    for(int i=0; i<samples; i++)
        ok = ok && approx(18, out[i], expect[i]);
    return ok;
}

TEST(minimal)
{
    // 2D 1x1 texture, no filtering
    texture tex(2, 1, texture::UBYTE);
    char data = 123;
    int wh[] = { 1, 1 };
    tex.set_mip(0, wh, &data);

    void* mems[] = { tex.get_mem() };
    const int memidx = 0;

    vecgen vg(optstring().c_str());
    vr texel, coords[2];
    coords[0] = vg.imm(0.5);
    coords[1] = vg.imm(0.5);
    tex.sample(&vg, memidx, coords, 0, 0, &texel);
    vg.output(0) = texel;

    float out;
    vecgen_run(vg, 0, 0, &out, mems);
    return approx(18, out, 123/255.0);
}

// Cube of three-component values that recapitulate the coordinates.
// No filtering.  Uses a border color and checks overflow.
TEST(unfiltered3d)
{
    int i, j, k;
    float out[3];
    const int sz = 13;

    texture tex(3, 4, texture::USHORT);
    tex.set_dims(false, false, texture::BORDER);
    float border[] = { -1, -1, -1 };
    tex.set_border_color(border);

    unsigned short *buf = new unsigned short[4*sz*sz*sz];
    for(i=0; i<sz; i++) for(j=0; j<sz; j++) for(k=0; k<sz; k++) {
        unsigned short *p = &buf[4*((k*sz+j)*sz+i)];
        p[0] = i; p[1] = j; p[2] = k; p[3] = 0xffff;
    }
    int extents[] = { sz, sz, sz };
    tex.set_mip(0, extents, (char*)buf);

    void *mems[] = { tex.get_mem() };
    const int memidx = 0;

    vecgen vg(optstring().c_str());
    vr coords[] = { vg.input(0), vg.input(1), vg.input(2) };
    vr outs[] = { vg.output(0), vg.output(1), vg.output(2), vg.scratch() };
    tex.sample(&vg, memidx, coords, 0, 0, outs);
    vg.codegen();

    // Spread 1 texel in each direction to check border behavior
    bool ok = false;
    for(i=-1; i<=sz; i++) for(j=-1; j<=sz; j++) for(k=-1; k<=sz; k++) {
        const float off = 0.5; // should manually try differnt values
        float in[] = { (i+off)/sz, (j+off)/sz, (k+off)/sz };
        vecgen_run_no_gen(vg, 0, in, out, mems);

        unsigned short *p = &buf[4*((k*sz+j)*sz+i)];
        bool border = i < 0 || i >= sz || j < 0 || j >= sz || k < 0 || k >= sz;
        for(int x=0; x<3; x++) {
            float expect = border ? -1 : p[x] * (1.0/0xffff);
            if(!approx(18, out[x], expect))
                goto done;
        }
    }
    ok = true;

 done:
    if(!ok)
        printf("ERR: i %d j %d k %d = (%f, %f, %f)\n", i, j, k, out[0], out[1], out[2]);

    delete[] buf;
    return ok;
}

// Emit code for a "typical" 2D trilinear 32 bit use case for hand
// optimization/analysis.  Always succeeds.
TEST(typical)
{
    texture tex(2, 4, texture::UBYTE);
    //tex.set_max_aniso_samples(16);
    tex.set_dims(true, false, texture::REPEAT);
    tex.set_mip_mode(texture::TRILINEAR);
    int mip0[16], mip1[4], mip2[1];
    int wh[2];
    wh[0] = wh[1] = 4;
    tex.set_mip(0, wh, (char*)mip0);
    wh[0] = wh[1] = 2;
    tex.set_mip(1, wh, (char*)mip1);
    wh[0] = wh[1] = 1;
    tex.set_mip(2, wh, (char*)mip2);

    void *mems[] = { tex.get_mem() };
    const int memidx = 0;

    vecgen vg(optstring().c_str());
    vr coords[] = { vg.input(0), vg.input(1), vg.input(2) };
    vr outs[] = { vg.output(0), vg.output(1), vg.output(2), vg.output(3) };
    vr ddx[] = { vg.imm(.25), vg.imm(0) };
    vr ddy[] = { vg.imm(0), vg.imm(.25) };
    tex.sample(&vg, memidx, coords, ddx, ddy, outs);

    float in[3], out[4];
    vecgen_run(vg, 0, in, out, mems);
    return true;
}
