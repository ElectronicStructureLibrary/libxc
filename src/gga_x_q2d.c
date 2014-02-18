/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

#define XC_GGA_X_Q2D          48 /* Chiodo et al  */

static void
gga_x_q2D_init(XC(func_type) *p)
{
  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_X_PBE_SOL, p->nspin);

}

void XC(gga_x_q2d_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static const FLOAT cc  = 100, c1 = 0.5217;

  FLOAT ss, ss2, ss4, ss6, ss_2, a, da, d2a, d3a;
  FLOAT f1, df1, d2f1, d3f1, f2, df2, d2f2, d3f2;

  ss   = X2S*x;
  ss2  = ss*ss;
  ss4  = ss2*ss2;
  ss6  = ss2*ss4;
  ss_2 = SQRT(ss);

  XC(gga_x_pbe_enhance)(p->func_aux[0], order, x, &a, &da, &d2a, &d3a);

  f1 = a*(cc - ss4) + c1*ss*ss2*ss_2*(1.0 + ss2);
  f2 = cc + ss6;
  *f = f1/f2;
  
  if(order < 1) return;

  da /= X2S;

  df1 = da*(cc - ss4) - 4.0*a*ss*ss2 + c1*ss2*ss_2*(7.0 + 11.0*ss2)/2.0;
  df2 = 6.0*ss*ss4;

  *dfdx  = X2S*DFRACTION(f1, df1, f2, df2);

  if(order < 2) return;

  d2a /= X2S*X2S;

  d2f1 = d2a*(cc - ss4) - 8.0*da*ss*ss2 - 12.0*a*ss2 + c1*ss*ss_2*(35.0 + 99.0*ss2)/4.0;
  d2f2 = 30.0*ss4;

  *d2fdx2 = X2S*X2S*D2FRACTION(f1, df1, d2f1, f2, df2, d2f2);

  if(order < 3) return;

  d3a /= X2S*X2S*X2S;

  d3f1 = d3a*(cc - ss4) - 12.0*ss*ss2*d2a - 36.0*da*ss2 - 24.0*a*ss + 21.0*c1*ss_2*(5.0 + 33.0*ss2)/8.0;
  d3f2 = 120.0*ss*ss2;
  
  *d3fdx3 = X2S*X2S*X2S*D3FRACTION(f1, df1, d2f1, d3f1, f2, df2, d2f2, d3f2);
}

#define func XC(gga_x_q2d_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_q2d) = {
  XC_GGA_X_Q2D,
  XC_EXCHANGE,
  "Chiodo et al",
  XC_FAMILY_GGA,
  "L Chiodo, LA Constantin, E Fabiano, and F Della Sala, Phys. Rev. Lett. 108, 126402 (2012)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-23, 0.0, 1e-32,
  gga_x_q2D_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};
