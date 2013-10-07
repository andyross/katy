// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#include <vector>
#include "vecgen.hpp"
#include "raster.hpp"

using namespace std;

// Some theory first.  Each attribute is defined to be linear in 3D
// (X/Y/W -- not Z, which is the scaled fragment depth and acts here
// as an attribute; W is the linear eye space depth coordinate), that
// is:
//
//   att = A*x  + B*y + C*w + D   for some A/B/C/D
//
// But D can be zero by construction.  Why?  The interpolation over
// the triangle is constrained only in the 2D plane of the triangle --
// the family of linear functions 3D functions that match it can be
// sheared outside the plane in any direction we choose.  So we simply
// choose the shear to draw a line-of-constant-attribute from the
// origin to the point on the plane where att is zero.  This works as
// long as the origin is not on the plane (but that's OK because such
// a case is degenerate: the triangle is edge-on).  Thus:
//
//   att = A*x + B*y * C*w
//
// Recognize that we're interpolating in post-perspective-divide
// screen space using "sx" proportional to x/w and "sy" that goes as
// y/w.  Rescaling A/B/C we can then write:
//
//     att = A*sx*w + B*sy*w + C*w
//         = w * (A*sx + B*sy + C)
//   att/w = A*sx + B*sy + C
//
// That is, "att/w" is linear in screen space.  We don't even need to
// solve for A/B/C: For each attribute-over-w "Q", the 2D attribute steps
// can be gotten via a 2x2 matrix inverse (which is dependent only on
// coordinates and can be computed once, not per-attribute):
//
//   (bq-aq) = (bx-ax)*dqdx + (by-ay)*dqdy
//   (cq-aq) = (cx-ax)*dqdx + (cy-ay)*dqdy
//
//   |bx-ax by-ay|   |dqdx|   |bq-aq|
//   |cx-ax cy-ay| * |dqdy| = |cq-aq|
//
//   |dqdx|     |bx-ax by-ay|-1   |bq-aq|
//   |dqdy|  == |cx-ax cy-ay|   * |cq-aq|


static void swapif(const vr &test, vr &a, vr &b)
{
    vr a2=a;
    a = mux(test, b, a);
    b = mux(test, a2, b);
}

void rasterize_tri(vecgen *vg, int viewx, int viewy, int vieww, int viewh,
                   vr *a, vr *b, vr *c, int natt_in, vr *aatt_in, vr *batt_in, vr *catt_in,
                   frag_fn fragment, void *frag_arg)
{
    // Might be immediates
    for(int i=0; i<4; i++) { a[i].setvg(vg); b[i].setvg(vg); c[i].setvg(vg);  }
    for(int i=0; i<natt_in; i++) { aatt_in[i].setvg(vg); batt_in[i].setvg(vg); catt_in[i].setvg(vg); }

    // Renaming
    vr ax=a[0], ay=a[1], az=a[2], aw=a[3];
    vr bx=b[0], by=b[1], bz=b[2], bw=b[3];
    vr cx=c[0], cy=c[1], cz=c[2], cw=c[3];

    // flip Y sign to convert cartesian to raster
    ay ^= imm_i(0x80000000);
    by ^= imm_i(0x80000000);
    cy ^= imm_i(0x80000000);

    // Perspective division.  Note we use a coordinate "m" to
    // represent 1/w (get it?!).
    vr am = recip2(aw); ax *= am; ay *= am; az *= am;
    vr bm = recip2(bw); bx *= bm; by *= bm; bz *= bm;
    vr cm = recip2(cw); cx *= cm; cy *= cm; cz *= cm;

    // Convert [-1:1] to [viewx/y,viewx/y+vieww/h] screen coordinates,
    // placing the origin at the center (!) of the bottom/left pixel.
    ax = ax * (vieww/2) + ((vieww/2.) + viewx - 0.5);
    ay = ay * (viewh/2) + ((viewh/2.) + viewy - 0.5);
    bx = bx * (vieww/2) + ((vieww/2.) + viewx - 0.5);
    by = by * (viewh/2) + ((viewh/2.) + viewy - 0.5);
    cx = cx * (vieww/2) + ((vieww/2.) + viewx - 0.5);
    cy = cy * (viewh/2) + ((viewh/2.) + viewy - 0.5);

    // Assemble a set of attributes (stored as att/w): include Z and
    // M, even though they interpolate somewhat differently, just to
    // have a place to put them.
    int natt = natt_in + 2, zi = natt-2, mi = zi+1;
    vector<vr> aatt, batt, catt;
    for(int i=0; i<natt_in; i++) {
        aatt.push_back(aatt_in[i] * am);
        batt.push_back(batt_in[i] * bm);
        catt.push_back(catt_in[i] * cm);
    }
    aatt.push_back(az); aatt.push_back(am);
    batt.push_back(bz); batt.push_back(bm);
    catt.push_back(cz); catt.push_back(cm);

    // Compute the matrix inverse for the 2D att/w derivatives
    vr m00 = cy-ay, m01 = ay-by;
    vr m10 = ax-cx, m11 = bx-ax;
    vr idet = recip2(m00*m11 - m01*m10);
    m00 *= idet; m01 *= idet; m10 *= idet; m11 *= idet;

    // Build the 2D derivatives
    vector<vr> ddx(natt), ddy(natt);
    for(int i=0; i<zi; i++) {
        vr ba = batt[i] - aatt[i], ca = catt[i] - aatt[i];
        ddx[i] = ba*m00 + ca*m01;
        ddy[i] = ba*m10 + ca*m11;
    }
    ddx[zi] = (bz-az)*m00 + (cz-az)*m01;
    ddy[zi] = (bz-az)*m10 + (cz-az)*m11;
    ddx[mi] = (bm-am)*m00 + (cm-am)*m01;
    ddy[mi] = (bm-am)*m10 + (cm-am)*m11;

    // Now we can interpolate anywhere on the plane given only the 2D
    // screen coordinates, so figure out the pixels to rasterize.
    // First sort the 2D points by Y into top/mid/bottom coordinates:
    vr topx=ax, midx=bx, botx=cx, topy=ay, midy=by, boty=cy;
    vr sw0 = topy > midy; swapif(sw0, topx, midx); swapif(sw0, topy, midy);
    vr sw1 = topy > boty; swapif(sw1, topx, botx); swapif(sw1, topy, boty);
    vr sw2 = midy > boty; swapif(sw2, midx, botx); swapif(sw2, midy, boty);

    // And a left/right edge line (expressed as x=q*y+o) for each
    // half-triangle.  The left edge is the one with a lower q slope.
    // NOTE: the slope NaNs out (recip2 doesn't do Inf properly) when
    // faced with a horizontal upper or lower edge, that will still
    // work as long as the code never "uses" the zero-size half
    // triangle (i.e. get the < vs. <= tests right).  FIXME: MUST TEST
    // THAT CONDITION.
    vr botlq = (midx-botx) * recip2(midy-boty);
    vr botlo = botx - botlq * boty;
    vr botrq = (topx-botx) * recip2(topy-boty);
    vr botro = botx - botrq * boty;
    vr toplq = (topx-midx) * recip2(topy-midy);
    vr toplo = topx - toplq * topy;
    vr toprq = botrq;
    vr topro = botro;
    // Sort left/right:
    vr swapedge = toprq < toplq;
    swapif(swapedge, botro, botlo); swapif(swapedge, botrq, botlq);
    swapif(swapedge, topro, toplo); swapif(swapedge, toprq, toplq);

    // Now loop over spans from the top down, we're doing 2x2 fragment
    // blocks at a time over {x,y}{0,1}, where 0,0 is the bottom left
    vr y0 = max(ceil(topy), viewy);
    vr maxy = min(boty, viewy+viewh);
    VRLOOP(*vg) {
        vr y1 = y0 + 1;

        // Figure out the q/o values for the left/right sides of each
        // span and compute span boundaries.
        vr ql0 = mux(y0 < midy, toplq, botlq), ol0 = mux(y0 < midy, toplo, botlo);
        vr qr0 = mux(y0 < midy, toprq, botrq), or0 = mux(y0 < midy, topro, botro);
        vr ql1 = mux(y1 < midy, toplq, botlq), ol1 = mux(y1 < midy, toplo, botlo);
        vr qr1 = mux(y1 < midy, toprq, botrq), or1 = mux(y1 < midy, topro, botro);
        vr xmin0 = ql0 * y0 + ol0, xmax0 = qr0 * y0 + or0;
        vr xmin1 = ql1 * y1 + ol1, xmax1 = qr1 * y1 + or1;

        // clamp to viewport.
        xmin0 = max(xmin0, viewx-.5);  xmax0 = min(xmax0, viewx+vieww-.5);
        xmin1 = max(xmin1, viewx-.5);  xmax1 = min(xmax1, viewx+vieww-.5);

        vr x0 = floor(min(xmin0, xmin1)), xmax = max(xmax0, xmax1);

        // Compute attributes for the bottom left corner
        vector<vr> att00(natt);
        for(int i=0; i<natt; i++)
            att00[i] = aatt[i] + (x0-ax)*ddx[i] + (y0-ay)*ddy[i];

        VRLOOP(*vg) {
            vr x1 = x0 + 1;

            // Figure the other three fragments' attributes:
            vector<vr> att01(natt), att10(natt), att11(natt);
            for(int i=0; i<natt; i++) {
                att01[i] = att00[i] + ddx[i]; // bottom right
                att10[i] = att00[i] + ddy[i]; // top left
                att11[i] = att01[i] + ddy[i]; // top right
            }

            // Compute W for each fragment and multiply by att/w
            // (except for Z of course) to get the true attribute
            // value.  Save off att01 pre-munging to use for computing
            // the next iteration.
            vector<vr> att01_copy = att01;
            vector<vr> *atts[] = { &att00, &att01, &att10, &att11 };
            for(int i=0; i<4; i++) {
                vr w = recip((*atts[i])[mi]); // per-pixel recip is the fast/approximate version
                for(int j=0; j<zi; j++)
                    (*atts[i])[j] *= w;
            }

            // Four fragments:
            for(int i=0; i<4; i++) {
                vector<vr> rec, *att = atts[i];
                int top = !(i & 2), right = !!(i & 1);

                // Build the attribute/derivatives array (3 entries
                // each: value/ddx/ddy)
                for(int j=0; j<mi; j++) {
                    vr jddx = (*atts[ 1 + (top<<1) ])[j] - (*atts[ 0 + (top<<1) ])[j];
                    vr jddy = (*atts[ 2 +   right  ])[j] - (*atts[ 0 +   right  ])[j];
                    rec.push_back((*att)[j]);
                    rec.push_back(jddx);
                    rec.push_back(jddy);
                }

                // Fragment inclusion check vs. viewport and Zfar.
                vr x = right ? x1 : x0, y = top ? y0 : y1, z = (*att)[zi];
                vr span0 = top ? xmin0 : xmin1, span1 = top ? xmax0 : xmax1;
                vr dofrag = ((x >= span0 && x <= span1) &&
                             (top ? z <= 1 : (z <= 1 && y < maxy)));
                VRIF(*vg, dofrag)
                    fragment(x, y, &rec[0], rec.size(), frag_arg);
            }

            // Update bottom left "seed" for next fragment quad
            for(int i=0; i<natt; i++)
                att00[i] = att01_copy[i] + ddx[i];

            x0 += 2;
            VRIF(*vg, x0 >= xmax)
                vg->break_loop();
        }

        y0 += 2;
        VRIF(*vg, y0 > maxy)
            vg->break_loop();
    }
}
