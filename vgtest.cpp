// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include "math.hpp"
#include "test.hpp"

TEST(test1)
{
    vecgen vg(optstring().c_str());

    vg.output(0) = vg.input(0) + 1;
    vg.output(0) *= vg.constant(0) + i2f(vg.imm_i(7));

    float in, out, c=100;
    in = 7;
    vecgen_run(vg, &c, &in, &out, 0);
    return out == 856;
}

TEST(test2)
{
    vecgen vg(optstring().c_str());

    vg.output(0) = 7;
    VRLOOP(vg) {
        VRIF(vg, vg.input(0)) {
            vg.break_loop();
        }
    }

    float in, out;
    in = 1;
    vecgen_run(vg, 0, &in, &out, 0);
    return out == 7;
}

TEST(ifelse)
{
    vecgen vg(optstring().c_str());

    VRIF(vg, vg.input(0)) {
        vg.output(0) = 123;
    } VRELSE(vg) {
        vg.output(0) = 321;
    }

    float in, out=0;
    in = 0;
    vecgen_run(vg, 0, &in, &out, 0);
    return out == 321;
}

// Computes IN0 * IN1 (for integers >=1, should probably fix it to
// handle zero...)
TEST(itermul)
{
    vecgen vg(optstring().c_str());

    vr i, j;
    i = vg.input(0);
    vg.output(0) = 0;
    VRLOOP(vg) {
        j = vg.input(1);
        VRLOOP(vg) {
            vg.output(0) += 1;
            j -= 1;
            if(log_execution) {
                vg.log("[i, j] = ", i, j);
                vg.log("[out0] = ", vg.output(0));
            }
            VRIF(vg, j <= 0)
                vg.break_loop();
        }
        i -= 1;
        VRIF(vg, i <= 0)
            vg.break_loop();
    }

    float in[] = { 6, 3 }, out;
    vecgen_run(vg, 0, in, &out, 0);
    return out == 18;
}

TEST(loadstore)
{
    vecgen vg(optstring().c_str());

    // Read the IN0'th byte and put it in OUT0
    vg.output(0) = load(U8, 0, vg.input(0));

    // Put an "X" in its place
    store(U8, 0, vg.input(0), vg.imm_i('X'));

    int in=3, out;
    char s[5];
    strcpy(s, "spot");
    void* memarray[1];
    memarray[0] = s;
    vecgen_run(vg, 0, (float*)&in, (float*)&out, memarray);
    return out == 't' && strcmp(s, "spoX") == 0;
}

TEST(memsz)
{
    vecgen vg(optstring().c_str());

    vr tmp = load(S16, 0, vg.imm_i(1));
    store(S16, 0, vg.imm_i(0), tmp);

    vr tmp2 = load(U32, 1, vg.imm_i(1));
    store(U32, 1, vg.imm_i(0), tmp2);

    short shorts[2] = { -1, 99 };
    int ints[2] = { -1, 99 };
    void* memarray[2] = { shorts, ints };
    vecgen_run(vg, 0, 0, 0, memarray);
    return shorts[0] == 99 && ints[0] == 99;
}

// Test automatic conversion between floats and AVX booleans
TEST(boolconv)
{
    vecgen vg(optstring().c_str());

    vr zero = vg.imm(0);
    vr one = vg.imm(1);
    vr two = vg.imm(2);
    vr three = vg.imm(3);

    vg.output(0) = (one < two) + (three > zero); // should == 2.0
    vg.output(1) = mux(zero, three, two); // should == 2.0

    float out[2];
    vecgen_run(vg, 0, 0, out, 0);
    return out[0] == 2 && out[1] == 2;
}

// Forces register spills and makes sure they work
TEST(regspill)
{
    const int n = 16;
    vecgen vg(optstring().c_str());
    vr one = vg.input(0);

    vector<vr> regs;
    regs.push_back(vg.imm(0));
    for(int i=0; i<n; i++) {
        regs.push_back(vg.scratch());
        regs.back() = regs[regs.size()-2] + one;
    }

    vg.output(0) = 0;
    int total=0;
    for(unsigned int i=0; i<regs.size(); i++) {
        vg.output(0) += regs[i];
        total += i;
    }

    float in=1, out;
    vecgen_run(vg, 0, &in, &out, 0);
    return out == total;
}

// Bland operator code coverage test.  Doesn't hit everything, just
// ones that gcov told me are uncovered.
//   RCP RSQ OR FEQ FNE FGE FLOOR CEIL F2I FMAX FMIN BITANDN BITAND
//   BITOR BITXOR FMA
TEST(ops)
{
    vecgen vg(optstring().c_str());
    map<int,float> expect;

#define O(n) vg.output(n)
#define I(v) vg.imm(v)
#define B(v) vg.imm_i(v)
    int n=0;
    O(n) = rsqrt(I(765.4321));     expect[n++] = 1/sqrt(765.4321);
    O(n) = recip(I(765.4321));     expect[n++] = 1/765.4321;
    O(n) = I(0) || I(1);           expect[n++] = 0 || 1;
    O(n) = I(0) || I(0);           expect[n++] = 0 || 0;
    O(n) = I(1) || I(1);           expect[n++] = 1 || 1;
    O(n) = I(1) || I(0);           expect[n++] = 1 || 0;
    O(n) = I(2) == I(3);           expect[n++] = 2 == 3;
    O(n) = I(2) == I(2);           expect[n++] = 2 == 2;
    O(n) = I(2) != I(3);           expect[n++] = 2 != 3;
    O(n) = I(2) != I(2);           expect[n++] = 2 != 2;
    O(n) = I(0) >= I(1);           expect[n++] = 0 >= 1;
    O(n) = I(1) >= I(1);           expect[n++] = 1 >= 1;
    O(n) = I(2) >= I(1);           expect[n++] = 2 >= 1;
    O(n) = floor(I(2.0001));       expect[n++] = floorf(2.0001);
    O(n) = ceil(I(2.0001));        expect[n++] = ceilf(2.0001);
    O(n) = max(I(4), I(2));        expect[n++] = max(4, 2);
    O(n) = min(I(4), I(2));        expect[n++] = min(4, 2);
    O(n) = fma(I(1), I(10), I(3)); expect[n++] = 1+10*3;

    union { float f; int i; } a, b, c;
    a.i = 0x80000000; // negative zero
    b.i = (127<<24); // 1.0
    O(n) = f2i(I(12.5));     c.i = (int)12.5;  expect[n++] = c.f; // 0x0000000c
    O(n) = B(a.i) | B(b.i);  c.i = a.i | b.i;  expect[n++] = c.f; // -1
    O(n) = B(a.i) & B(b.i);  c.i = a.i & b.i;  expect[n++] = c.f; // +0
    O(n) = B(a.i) ^ B(b.i);  c.i = a.i ^ b.i;  expect[n++] = c.f; // -0

    O(n) = bitandnot(B(a.i),B(b.i)); c.i = a.i & ~b.i;  expect[n++] = c.f; // -0

    float *out = new float[n];
    vecgen_run(vg, 0, 0, out, 0);

    bool ok = true;

    // First two are approximations, check for ~11 bits of precision,
    // which is what Intel specifies for RCPPS/RSQRPS
    if(1/fabs(1 - out[0]/expect[0]) < 2000) ok = false;
    if(1/fabs(1 - out[1]/expect[1]) < 2000) ok = false;

    for(int i=2; i<n; i++)
        if(out[i] != expect[i])
            ok = false;
    delete[] out;
    return ok;
}

TEST(loopnest)
{
    vecgen vg(optstring().c_str());
    VRLOOP(vg) {
        VRLOOP(vg) {
            vg.output(0) = 19;
            vg.break_loop();
        }
        vg.break_loop();
    }

    float out;
    vecgen_run(vg, 0, 0, &out, 0);
    return out == 19;
}

TEST(loopcount)
{
    vecgen vg(optstring().c_str());

    // Count in0 down to 4.
    vg.output(0) = vg.input(0);
    VRLOOP(vg) {
        vg.output(0) -= 1;
        VRIF(vg, vg.output(0) == 4)
            vg.break_loop();
    }

    float in = 16, out;
    vecgen_run(vg, 0, &in, &out, 0);
    return out == 4;
}

TEST(loophoist)
{
    vecgen vg(optstring().c_str());
    VRLOOP(vg) {
        vg.break_loop(); // empty loop
    }
    VRLOOP(vg) {
        VRIF(vg, vg.input(0) > vg.input(1)) {
            vg.output(0) = 14;
            vg.break_loop();
        }
    }

    float out=12, in[2] = { 100, 90 };
    vecgen_run(vg, 0, in, &out, 0);
    return out == 14;
}

TEST(mux)
{
    vecgen vg(optstring().c_str());

    // OUT0 = IN0 ? IN1 : IN2
    vg.output(0) = mux(vg.input(0), vg.input(1), vg.input(2));

    // OUT1 = IN3 ? IN1 : IN2
    vg.output(1) = mux(vg.input(3), vg.input(1), vg.input(2));

    float in[4], out[2];
    in[0] = 0;
    in[1] = 12;
    in[2] = -99;
    in[3] = 1;
    vecgen_run(vg, 0, in, out, 0);
    return ((out[0] == ((in[0] != 0) ? in[1] : in[2])) &&
            (out[1] == ((in[3] != 0) ? in[1] : in[2])));
}

// Exercise a bunch of known compiler errors for code coverage
TEST(compileerrors)
{
#define CHECK(msg) catch(domain_error e) { \
        died = true;                       \
        if(strcmp(e.what(), (msg)))        \
            return false;                  \
    } if(!died) return false;

    bool died = false;
    try {
        vecgen vg(optstring().c_str());
        vg.break_loop();
    }
    CHECK("break outside of loop");

    try {
        vecgen vg(optstring().c_str());
        vg.end_loop();
    } CHECK("end_loop outside of loop");

    try {
        vecgen vg(optstring().c_str());
        vg.start_else();
    } CHECK("start_else outside of if");

    try {
        vecgen vg(optstring().c_str());
        vg.end_if();
    } CHECK("end_if outside of if");

    try {
        vecgen vg(optstring().c_str());
        vg.start_loop();
        vg.codegen();
    } CHECK("codegen with unclosed blocks");

    return true;
}

static int cullcb_iter;
void cull_bound(void* arg)
{
    (void)arg;
    cullcb_iter++;
}

TEST(cull)
{
    const int N = 25;
    vecgen vg(optstring().c_str());

    vr a, b, c;
    a = vg.imm(0);
    b = vg.imm(1);
    c = vg.imm(2);
    for(int i=0; i<N; i++) {
        a += 1;
        b += 1;
        c += 1;
        vr crec[] = { a, b, c };
        int id = -1;
        VRIF(vg, (i & 1) ? 0 : 1)
            id = vg.cull(crec, 3, 0, 1);
        vg.set_cull_bound(id, 3, 0, 2);
    }

    float out[3*(8*((N+7)/8))];
    void *mems[] = { out, (void*)0, reinterpret_cast<void*>(cull_bound), out+(3*8) };
    cullcb_iter = 0;
    vecgen_run(vg, 0, 0, 0, mems);

    long ocnt = (long)mems[1];
    int nrecs = ((float*)mems[0] - out)/3 + ocnt;
    if(nrecs != (N+1)/2)
        return false;

    if(cullcb_iter != nrecs-4)
        return false;

    {
        float a=0, b=1, c=2, *orec = out;
        int n = 0;
        for(int i=0; i<N; i++) {
            a++; b++; c++;
            if((i & 1) == 0) {
                if(orec[0+n] != a || orec[8+n] != b || orec[16+n] != c)
                    return false;
                if(++n == 8) {
                    orec += 3*8;
                    n = 0;
                }
            }
        }
    }

    return true;
}

TEST(count) {
    vecgen vg(optstring().c_str());

    const int count_idx = 0, instride = 1, outstride = 0;
    vg.set_count(count_idx, instride, outstride);

    // Generate an incrementing count, store it in the unstrided
    // output (note that this exercises the fact that "output"
    // registers have initial values and can be read).  Check the
    // count vs. a pre-cooked input array to be sure the strides are
    // working.
    VRIF(vg, vg.input(0) == vg.output(0))
        vg.output(0) += 1;

    const int iterations = 10;
    float *in = new float[iterations];
    for(int i=0; i<iterations; i++) {
        // Note we have to hack creation of all the "inputs" so
        // vecgen_run() sees a big array (one stride per iteration)
        // instead of just a single register.
        (void)vg.input(i);
        in[i] = i;
    }

    float out[] = { 0 };
    void* mem[] = { (void*)iterations };

    vecgen_run(vg, 0, in, out, mem);
    delete[] in;

    // Note that the interpreter doesn't support counted iteration and
    // will have only executed the first loop.  Pass anyway.
    return !avxmode || (mem[0] == 0 && out[0] == iterations);
}

TEST(predication) {
    vecgen vg(optstring().c_str());

    VRIF(vg, vg.input(0) != vg.output(0))
        vg.output(0) += 1;

    float in[] = { 0 };
    float out[] = { 0 };
    vecgen_run(vg, 0, in, out, 0);
    return out[0] == 0;
}

static bool approx(int bits, double a, double b) {
    const double prec = 1./(1<<bits);
    return (a == 0 && b == 0) || (b != 0 && fabs(1 - a/b) < prec);
}

TEST(math) {
    vecgen vg(optstring().c_str());

    vr f;
    f = vg.imm(1.3);
    vg.output(0) = log2(f);
    vg.output(1) = exp2(f);
    vr out2 = vg.output(2), out3 = vg.output(3);
    sincos(f, out2, out3);
    vg.output(4) = atan(f);

    float out[5];
    vecgen_run(vg, 0, 0, out, 0);
    bool ok = (approx(23, out[0], log2(1.3)) ||
               approx(23, out[1], exp2(1.3)) ||
               approx(23, out[2], sin(1.3)) ||
               approx(23, out[3], cos(1.3)) ||
               approx(23, out[4], atan(1.3)));
    return ok;
}

TEST(atan) {
    // Regression test for a AVX-only discontinuity in atan() due to a
    // misfeature in rcpps (will return exactly zero for values >=
    // 2^126)
    vecgen vg(optstring().c_str());
    vr x; x = vg.input(0);
    vg.output(0) = atan(x);
    vg.codegen();

    union { int i; float f; } in, out;
    bool ok = true;

    in.i = 0x7f7ffffe;
    vecgen_run_no_gen(vg, 0, &in.f, &out.f, 0);
    ok = ok && approx(22, out.f, 1.5707963);
    //printf("atan(%x/%g) = %x/%g\n", in.i, in.f, out.i, out.f);

    in.i++;
    vecgen_run_no_gen(vg, 0, &in.f, &out.f, 0);
    ok = ok && approx(22, out.f, 1.5707963);
    //printf("atan(%x/%g) = %x/%g\n", in.i, in.f, out.i, out.f);

    return ok;
}
