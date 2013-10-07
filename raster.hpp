// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _RASTER_HPP
#define _RASTER_HPP

class vecgen;
class vr;

// Triangle rasterization utility.  It's a function that generates
// into a vecgen, providing a per-fragment callback.

// Handles one rasterized fragment out of the inner loop.  sx/sy are
// the screen-space pixel coordinates (floating point, but with
// guaranteed-integer values).  The attribute array is in groups of
// three: the value of each attribute, followed by ddx then ddy.  The
// final attribute in the array is the Z coordinate in the range
// [0:1].
typedef void (*frag_fn)(vr sx, vr sy, vr *attdxdy, int natt, void *arg);

// a/b/c are the x/y/z/w coordinates of the input points.  a/b/catt
// are the arrays of attributes.  The framebuffer has a width and
// height, and potentially separate strides (in pixels, not bytes) for
// the color and depth buffers.  Framebuffer is assumed to be
// raster-ordered, not cartesian (i.e. y=0 is the top).  Likewise the
// viewport is in raster/y-down coordinates (i.e. *not* like
// glViewport()!) and assumed to be valid (i.e. entirely contained
// within the framebuffer).
void rasterize_tri(vecgen *vg, int viewx, int viewy, int vieww, int viewh,
                   vr *a, vr *b, vr *c, int natt_in, vr *aatt_in, vr *batt_in, vr *catt_in,
                   frag_fn fragment, void *frag_arg);

#endif // _RASTER_HPP
