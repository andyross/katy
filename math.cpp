// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include "math.hpp"

// The implementations here are based on the SunPro code used
// basically everywhere, generally with reduced polynomials or
// simplified functions that work better in float precision.

// Other stuff to consider adding:
//
// cosh/sinh/acosh/asinh
// erf/tgamma/lgamma/j0/j1/jn
// expm1
// isfinite/isinf/isnan/logb/significand
// ldexp/scalbn
// fmod/remainder
// cbrt/hypot
// rint (do in vecgen)

static float pow2(int p)
{
    union { int i; float f; } u;
    u.i = (p+127) << 23;
    return u.f;
}

// Computes sin() in the range -pi/4:pi/4
static vr sin_kernel(vr x)
{
    vr S1  = -1.6666667e-1f;
    vr S2  =  8.3333333e-3f;
    vr S3  = -1.9841270e-4f;
    vr s = x & imm_i(0x80000000); // store original sign
    vr x2 = x * x;
    x &= imm_i(0x7fffffff);       // mask off sign
    return s | (x + (x * x2*(S1 + x2*(S2 + x2*S3))));
    return s | (x + (x * x2*(S1 + x2*(S2 + x2*S3))));
}

// Computes cos() in the range -pi/4:pi/4
static vr cos_kernel(vr x)
{
    vr C1 =  4.1666667e-2f;
    vr C2 = -1.3888889e-3f;
    vr C3 =  2.4801587e-5f;
    x &= imm_i(0x7fffffff); // mask off sign
    vr x2  = x * x;
    return imm(1) - x2*(imm(0.5) - x2*(C1 + x2*(C2 + x2*C3)));
}

// Tested to 8 ulp vs. the input scale
void sincos(vr v, vr &sin_out, vr &cos_out)
{
    vr PI = 3.1415927;
    vr iTWOPI = 0.15915494; // 1/(2*PI)
    vr TWOPI = 6.2831853;   // 2*PI
    vr PIQ = 0.78539816;    // PI/4
    vr PIH = 1.5707963;     // PI/2
    vr PI3Q = 2.3561945;    // 3*PI/4

    // Fold input to positive numbers via cos(-x)=cos(x) and sin(-x) =
    // -sin(x) (i.e. store the sign bit of the input to mask into the
    // sine output).
    vr sin_sign = v & imm_i(0x80000000);
    v &= imm_i(0x7fffffff);

    // Reduce/modulus to [0:2pi]
    vr reduce = v < 0 || v > TWOPI;
    v -= mux(reduce, floor(v * iTWOPI) * TWOPI, 0);

    // Fold to [0:PI] by exploiting more symmetry: sin/cos(x-pi)=-sin/cos(x)
    vr fold = v > PI;
    v = mux(fold, v-PI, v);
    vr cos_sign = mux(fold, imm_i(0x80000000), imm_i(0)); // FIXME: noop in AVX
    sin_sign ^= cos_sign;

    // We have three ranges to test:
    // + "lo": If x < pi/4 we can compute directly
    // + "mid": If x is in [pi/4, 3*pi/4]:
    //    + x = x - pi/2
    //    + swap the sin/cos output
    //    + invert the cos sign
    // + "hi": If x > 3*pi/4
    //    + x = x - pi
    //    + invert both sign bits
    vr lo = v < PIQ;
    vr hi = v > PI3Q;
    v -= mux(lo, 0, mux(hi, PI, PIH));
    sin_sign ^= mux(hi, imm_i(0x80000000), imm_i(0)); // FIXME: noop in AVX
    cos_sign ^= mux(lo, imm_i(0), imm_i(0x80000000));

    vr s = sin_kernel(v);
    vr c = cos_kernel(v);
    sin_out = sin_sign ^ mux((lo || hi), s, c);
    cos_out = cos_sign ^ mux((lo || hi), c, s);
}

vr log2(vr x)
{
    const float TWO_OVER_LN2 = 2.885390f; // == 2/ln(2)
    vr ILN2_2 = TWO_OVER_LN2;
    vr SQRT2 = 1.414213f;
    vr C1 = TWO_OVER_LN2/3.0;
    vr C2 = TWO_OVER_LN2/5.0;
    vr C3 = TWO_OVER_LN2/7.0;

    // Split exponent & mantissa
    vr exp = i2f(x & imm_i(0x7f800000)) * pow2(-23) - 127;
    x &= imm_i(0x007fffff);
    x |= imm_i(127 << 23);

    // Reduce x from [1,2] to [0:sqrt(2)], bump exponent if needed
    vr shrink = (x > SQRT2);
    x   = mux(shrink, x*0.5,   x);
    exp = mux(shrink, exp+1, exp);

    // Compute polynomial approximation of log2(2**exp + x)
    vr s = (x-1) * recip2(x+1);
    vr s2 = s * s;
    return exp + s * (ILN2_2 + ((C3*s2 + C2)*s2 + C1)*s2);
}

vr exp2(vr x)
{
    vr LN2 = 0.69314718;
    vr C1 =  1/6.;
    vr C2 = -1/360.;
    vr C3 =  1/15120.;

    // Premultiply by ln(2) because we're computing a base2 exponent
    // (save the original as it avoids a division later though)
    vr x_over_ln2 = x;
    x *= LN2;

    // Compute k (an integer) and r such that x == k*ln(2) + r
    // So r is constrained to [0,ln(2)]
    vr k = floor(x_over_ln2);
    vr r = x - k * LN2;

    // Now exp(x) == exp(ln(2)*k)*exp(r) == pow2(k)*exp(r)
    // and we can compute exp(r) in this range via polynomial:
    vr pow2k = f2i((k + 127) * pow2(23));
    vr r2 = r * r;
    vr R = ((C3*r2 + C2)*r2 + C1)*r2 + 2;
    return pow2k * ((r*2)*recip2(R-r) + 1);
}

// Good to 3ULP vs. glibc 2.14.90
// Does not handle +/- Inf inputs (should map to +/-PI)
vr atan(vr x)
{
    vr C0 =  3.333333e-01;
    vr C1 = -2.000000e-01;
    vr C2 =  1.428571e-01;
    vr C3 = -1.111111e-01;
    vr C4 =  9.090887e-02;
    vr C5 = -7.691876e-02;
    vr C6 =  6.661073e-02;
    vr C7 = -5.833570e-02;
    vr C8 =  4.976878e-02;
    vr C9 = -3.653157e-02;
    vr C10 = 1.628582e-02;

    vr sign = x & imm_i(0x80000000);
    x &= imm_i(0x7fffffff);

    // This implementation eliminates most of the subranges in the
    // SunPro code because they don't add meaningful precision for a
    // float.
    //
    // Note the recip2 expression which is essentially computing
    // "1/(x*1.5+1)": there's an overflow issue on AVX here.  The
    // rcpps instruction will (not-quite-incorrectly, but very
    // inconveniently) return exactly zero for values with an exponent
    // of 2^126 or 2^127 (presumably because those results are not
    // representable without resorting to a denormalized value the
    // hardware doesn't support), which causes x to become identically
    // zero and breaks the polynomial which would otherwise work for
    // any finite value.  The multiplication by the extra constant of
    // 1/8 isn't needless: it's there to ensure that the argument to
    // recip is always less than 2^126, and thus the result is finite
    // for all finite inputs.
    vr lo = x < 0.52;
    vr off = mux(lo, 0, .9827937);
    x = mux(lo, x, ((x-1.5) * imm(0.125)) * recip2(x*.1875+0.125));

    vr x2 = x*x;
    vr x4 = x2*x2;
    vr s1 = x2*(C0 + x4*(C2 + x4*(C4 + x4*(C6 + x4*(C8 + x4*C10)))));
    vr s2 = x4*(C1 + x4*(C3 + x4*(C5 + x4*(C7 + x4*C9))));
    return sign | (off + x - x*(s1+s2));
}

