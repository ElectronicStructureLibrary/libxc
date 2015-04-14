/*
 Copyright (C) 2008 Georg Madsen

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_EV93  35 /* Engel and Vosko */

void XC(gga_x_ev93_enhance)
     (const XC(func_type) *p, int order, FLOAT x, 
      FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static FLOAT 
    a1 = 1.647127,
    a2 = 0.980118, 
    a3 = 0.017399, 
    b1 = 1.523671, 
    b2 = 0.367229, 
    b3 = 0.011282;

  FLOAT ss, ss2, ss4, ss6;
  FLOAT   den,   num;
  FLOAT  dden,  dnum;
  FLOAT d2den, d2num;
  FLOAT d3den, d3num;

  ss  = X2S*x;
  ss2 = ss*ss;
  ss4 = ss2*ss2;
  ss6 = ss4*ss2;

  num = 1.0 + a1*ss2 + a2*ss4 + a3*ss6;
  den = 1.0 + b1*ss2 + b2*ss4 + b3*ss6;

  *f = num/den;

  if(order < 1) return;

  dnum = ss*(2.0*a1 + 4.0*a2*ss2 + 6.0*a3*ss4);
  dden = ss*(2.0*b1 + 4.0*b2*ss2 + 6.0*b3*ss4);

  *dfdx  = DFRACTION(num, dnum, den, dden);
  *dfdx *= X2S;

  if(order < 2) return;

  d2num = 2.0*a1 + 4.0*3.0*a2*ss2 + 6.0*5.0*a3*ss4;
  d2den = 2.0*b1 + 4.0*3.0*b2*ss2 + 6.0*5.0*b3*ss4;
  
  *d2fdx2  = D2FRACTION(num, dnum, d2num, den, dden, d2den);
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3num = 4.0*3.0*2.0*a2*ss + 6.0*5.0*4.0*a3*ss2*ss;
  d3den = 4.0*3.0*2.0*b2*ss + 6.0*5.0*4.0*b3*ss2*ss;

  *d3fdx3  = D3FRACTION(num, dnum, d2num, d3num, den, dden, d2den, d3den);
  *d3fdx3 *= X2S*X2S*X2S;
}

#define func XC(gga_x_ev93_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_ev93) = {
  XC_GGA_X_EV93,
  XC_EXCHANGE,
  "Engel and Vosko",
  XC_FAMILY_GGA,
  {&xc_ref_Engel1993_13164, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL,
  NULL, NULL,
  work_gga_x,
  NULL
};
