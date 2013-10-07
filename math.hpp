// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _MATH_HPP
#define _MATH_HPP

#include "vecgen.hpp"

// Simple math library for vr's

// Basic functions:
vr atan(vr x);
vr exp2(vr x);
vr log2(vr x);
void sincos(vr v, vr &sin_out, vr &cos_out);

// Derived from the above:

inline vr atan2(vr x, vr y)
{
    const float PI = 3.1415927;
    const float TWOPI = 6.2831853;
    vr out = atan(y * recip2(x));
    out = mux(x < 0, out + PI, mux(y < 0, out + TWOPI, out));
    out = mux(x == 0, (y & imm_i(0x80000000)) | PI, out);
    return out;
}

inline vr sin(vr x) { vr a, b; sincos(x, a, b); return a; }
inline vr cos(vr x) { vr a, b; sincos(x, a, b); return b; }
inline vr tan(vr x) { vr a, b; sincos(x, a, b); return a * recip2(b); }
inline vr asin(vr x) { return atan(x * rsqrt2(imm(1)-x*x)); }
inline vr acos(vr x) { return atan(recip2(x * rsqrt2(imm(1)-x*x))); }
inline vr log(vr x) { return log2(x) * 0.69314718; }
inline vr exp(vr x) { return exp2(x * 1.4426950); }
inline vr log10(vr x) { return log2(x) * 2.3025851; }
inline vr exp10(vr x) { return exp2(x * 0.43429448); }
inline vr pow(vr a, vr b) { return exp2(b * log2(a)); }

#endif // _MATH_HPP
