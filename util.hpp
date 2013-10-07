// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _UTIL_HPP
#define _UTIL_HPP

typedef union { int i; float f; } uif;

// Returns a value that is less than the positive-definite (!) input
// by the smallest representable ("unit in last place") fraction.  The
// simple subtraction works because if it rolls over the mantissa, it
// borrows one from the exponent, which is what we want.  It even
// correctly rolls into denormalied numbers and hits positive zero
// before rolling into NaN space.
inline float ulp_less(double val)
{
    uif u;
    u.f = val;
    u.i -= 1;
    return u.f;
}

#endif // _UTIL_HPP
