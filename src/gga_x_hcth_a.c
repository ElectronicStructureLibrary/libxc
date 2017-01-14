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

#define XC_GGA_X_HCTH_A          34 /* HCTH-A */

static void 
func(const XC(func_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  const FLOAT beta = 0.0042, gamma = 6;
  const FLOAT c0 = 1.09878, c1 = -2.51173, c2 = 0.0156233;
  FLOAT x2, aux1, aux2, f1, f2, f3, df1, df2, df3, d2f1, d2f2, d2f3, d3f1, d3f2, d3f3;

  x2 = x*x;

  f1 = beta/X_FACTOR_C*x2;
  f2 = 1.0 + gamma*beta*x*ASINH(x);
  f3 = beta*f2*f2;
  *f = c0 + c1*f1/f2 + c2*f1/f3;

  if(order < 1) return;

  aux1 = 1.0 + x2;
  aux2 = SQRT(aux1);

  df1 = 2.0*beta/X_FACTOR_C*x;
  df2 = gamma*beta*(ASINH(x) + x/aux2);
  df3 = 2.0*beta*f2*df2;

  *dfdx = c1*DFRACTION(f1, df1, f2, df2) 
    + c2*DFRACTION(f1, df1, f3, df3);

  if(order < 2) return;

  d2f1 = 2.0*beta/X_FACTOR_C;
  d2f2 = gamma*beta*(2.0 + x2)/(aux1*aux2);
  d2f3 = 2.0*beta*(df2*df2 + f2*d2f2);

  *d2fdx2 = c1*D2FRACTION(f1, df1, d2f1, f2, df2, d2f2) 
    + c2*D2FRACTION(f1, df1, d2f1, f3, df3, d2f3);

  if(order < 3) return;

  d3f1 = 0.0;
  d3f2 = -beta*gamma*x*(4.0 + x2)/(aux1*aux1*aux2);
  d3f3 = 2.0*beta*(3.0*df2*d2f2 + f2*d3f2);

  *d3fdx3 = c1*D3FRACTION(f1, df1, d2f1, d3f1, f2, df2, d2f2, d3f2) 
    + c2*D3FRACTION(f1, df1, d2f1, d3f1, f3, df3, d2f3, d3f3);
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_hcth_a) = {
  XC_GGA_X_HCTH_A,
  XC_EXCHANGE,
  "HCTH-A",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, 
  NULL,
  NULL,
  work_gga_x,
  NULL
};
