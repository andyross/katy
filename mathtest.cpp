// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <sys/mman.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <list>
#include "thread.hpp"
#include "math.hpp"

// Uses vector code in parallel to compute the math functions across
// their full domains and check for accuracy vs. the native C library.
//
// Obviously spends most of its time in the scalar libm code, but
// still acts as a pretty good stress on the AVX layer.

using namespace std;

struct work_item { thread_fn fn; void * arg; sem *done; };
list<work_item> work_queue;
mutex work_lock;
sem work_sem, tests_sem;
int nworkers=0;

const int nrecs = 254; // number of 8-entry records: comes out to a 16k buffer

struct test {
    const char* name;
    vr (*gen_fn)(vr);
    double (*ref_fn)(double);
    float min, max;

    // The sin/cos/atan functions have zero crossings and a limited
    // range (e.g. +/-1).  A tiny difference in the x intercept can
    // produce exploding "error" terms, so for these we clamp the
    // denominator of the error term to the range instead of the
    // maximum absolute value of the two results.  Cheating?  Who's to
    // say.
    float range;
};
extern struct test tests[];

struct vec_item {
    vecgen_fn code;
    float *buf;
};

// Maximum representable positive quantity
float float_max()
{
    union { int i; float f; } u;
    u.i = 0x7f7fffff;
    return u.f;
}

void worker_thread(void *unused)
{
    (void)unused;
    work_lock.lock(); nworkers++; work_lock.unlock();
    while(1) {
        work_item wi;
        work_sem.down();
        work_lock.lock();
        wi = work_queue.front();
        work_queue.pop_front();
        work_lock.unlock();
        wi.fn(wi.arg);
        wi.done->up();
    }
}

void vec_work(vec_item *vi)
{
    float *outs = vi->buf, *ins = &vi->buf[4*8];
    void *mems[] = { (void*)nrecs };
    int mskmem[16]; // need 8, extra for alignment
    int *msk = (int*)(((long)&mskmem[8]) & (~31L));
    for(int i=0; i<8; i++)
        msk[i] = 0x80000000;

    vi->code(0, ins, outs, mems, (float*)msk);
    delete vi;
}

// A worker thread will do the item and up the semaphore when complete
void add_work_item(thread_fn fn, void* arg, sem* s)
{
    work_item wi = { fn, arg, s };
    work_lock.lock();
    work_queue.push_back(wi);
    work_lock.unlock();
    work_sem.up();
}

// Returns the "next highest" IEEE 32 bit normalized finite (or zero) float.
float ieee_next(float x)
{
    union { int i; float f; } u;
    u.f = x;
    if(u.i & 0x80000000) { // negative?
        if(u.i == (int)0x80000000) {
            u.i = 0; // negative zero -> positive zero
        } else {
            if((--u.i & 0x7f800000) == 0) // decrement
                u.i = 0x80000000; // rollover -> negative zero
        }
    } else if(u.i == 0) {
        u.i = 0x00800000; // positive zero -> smallest positive number
    } else {
        u.i++;
    }
    return u.f;
}

// Returns the scalar index of a single float in an AVX array of
// records nfields large, in the specified field (i.e. OUT3==3), at
// the specified index (where 0-7 are in the first AVX record, 8-15 in
// the second, etc...)
inline int fieldloc(int nfields, int field, int idx)
{
    int rec = idx >> 3, i = idx % 8;
    return 8*(nfields*rec + field) + i;
}

vr abs(vr a) { return a & imm_i(0x7fffffff); }

int to_int(float f)
{
    union { int i; float f; } u;
    u.f = f;
    return u.i;
}

void test_thread(test *t)
{
    // The vecgen works on a rolling array of input records (0: input
    // value, 1: libc expectation value) and a static output set
    // containing the maximum error for that thread (0: input, 1:
    // libc, 2: result, 3: maxerr).
    vecgen vg;
    vr x = vg.input(0);
    vr libc = vg.input(1);
    vr result = t->gen_fn(x);
    vr scale = max(abs(result), abs(libc));
    if(t->range)
        scale = max(t->range, scale);
    vr err; err = abs(libc - result) * recip2(scale);
    VRIF(vg, err > vg.output(3)) {
        vg.output(0) = x;
        vg.output(1) = libc;
        vg.output(2) = result;
        vg.output(3) = err;
    }
    vg.set_count(0, 2, 0); // memidx, in stride, out stride
    vg.codegen();

    void* code = mmap(0, vg.code_size(), PROT_EXEC|PROT_WRITE,
                      MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
    vg.get_code((char*)code);

    const int nbufs = 4 * nworkers;
    const int wbsz = 8* (2 * nrecs + 4);
    typedef float workbuf[wbsz];
    workbuf *bufs;
    posix_memalign((void**)&bufs, 32, 4*nbufs*wbsz);

    sem done_sem;
    int currbuf=0, n=0;
    float max_err=0, max_input=0, max_libc=0, max_result=0;
    for(float x = t->min; /* see below */; x = ieee_next(x)) {
        float *recs = &bufs[currbuf][8*4]; // skip the first four entries, those are output
        recs[fieldloc(2, 0, n)] = x;
        recs[fieldloc(2, 1, n)] = t->ref_fn(x);
        n++;
        if(n < 8*nrecs && x != t->max)
            continue;

        // Buffer full.  Submit.
        for(int i=0; i<4*8; i++)
            bufs[currbuf][i] = -1;

        vec_item vi = { (vecgen_fn)code, (float*)&bufs[currbuf] };
        add_work_item((thread_fn)vec_work, new vec_item(vi), &done_sem);
        currbuf++;
        n = 0;

        // All buffers queueud (or at end of data)?
        // FIXME: would work better if we had a separate "return
        // queue" we could read from instead of doing it in batches,
        // but that complicates the interface
        if(currbuf == nbufs || x == t->max) {
            for(int i=0; i<currbuf; i++)
                done_sem.down(); // Wait for completion

            // Check for maximum
            for(int i=0; i<nbufs; i++) {
                for(int j=0; j<8; j++) {
                    float err = bufs[i][fieldloc(0, 3, j)];
                    if(err > max_err) {
                        max_err = err;
                        max_input = bufs[i][fieldloc(0, 0, j)]; // input
                        max_libc = bufs[i][fieldloc(0, 1, j)]; // libc
                        max_result = bufs[i][fieldloc(0, 2, j)]; // result;
                    }
                }
            }
            currbuf = 0;
        }

        // Can't exit using < in the for() test because float_next(max) might be NaN
        if(x == t->max)
            break;
    }
    printf("%4s: max err %g (~%.1f ulp) at %g/0x%8.8x\n"
           "      (expect %g/0x%8.8x, got %g/0x%8.8x)\n",
           t->name, max_err, max_err*(1<<24), max_input,
           to_int(max_input), max_libc, to_int(max_libc),
           max_result, to_int(max_result));
    free(bufs);
    tests_sem.up();
}

int count_cpus()
{
    int ncpus = 0;
    ifstream cpuinfo("/proc/cpuinfo");
    while(cpuinfo.good()) {
        string l;
        getline(cpuinfo, l);
        ncpus += (l.find("processor") == 0);
    }
    return ncpus;
}

vr log2_gen(vr x) { return log2(x); }
vr exp2_gen(vr x) { return exp2(x); }
vr atan_gen(vr x) { return atan(x); }
vr sin_gen(vr x) { return sin(x); }
vr cos_gen(vr x) { return cos(x); }

struct test tests[] = {
    // log2()'s minimum must bump twice from zero, because 0x00800001
    // produces -Inf on glibc (not quite wrong, but it blows up our
    // error computation)
    { "log2", log2_gen, log2, ieee_next(ieee_next(0)), float_max() },

    // The implementation of exp2 has a similar cyclic precision loss
    // to sin/cos (though not as bad) where the inherent error is of
    // the same order as the input value.
    { "exp2", exp2_gen, exp2, -126, 126 },

    { "atan", atan_gen, atan, -float_max(), float_max(), 1.5707963 },

    // Because the rounding is done in float precision, sin/cos lose a
    // ulp every cycle beyond zero (libm does too, but it's in double
    // precision and thus the error is below detection for the first
    // 2^29 cycles).  Just test two.
    { "sin", sin_gen, sin, -6.3, 6.3, 1 },
    { "cos", cos_gen, cos, -6.3, 6.3, 1 },
};

int main()
{
    const int ncpus = count_cpus();
    for(int i=0; i<ncpus; i++)
        start_thread(worker_thread, 0);

    int ntests = sizeof(tests)/sizeof(tests[0]);
    for(int i=0; i<ntests; i++)
        start_thread((thread_fn)test_thread, &tests[i]);

    for(int i=0; i<ntests; i++)
        tests_sem.down();
    return 0;
}
