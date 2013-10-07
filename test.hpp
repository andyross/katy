// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _TEST_HPP
#define _TEST_HPP

#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "vecgen.hpp"

using namespace std;

struct test_rec {
    const char *name; bool (*fn)();
    test_rec(const char *n, bool(*f)());
};
#define TEST(N) bool(N)(); const static test_rec t_##N(#N, N); bool(N)()

string optstring();
void vecgen_run(vecgen &vg, float *consts, float *in, float *out, void **mems);
void vecgen_run_no_gen(vecgen &vg, float *consts, float *in, float *out, void **mems);

extern bool log_execution;
extern bool avxmode;

#endif // _TEST_HPP
