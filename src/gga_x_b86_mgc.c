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
#include <assert.h>
#include "util.h"

#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalfa,beta,gamma (with mod. grad. correction) */

void XC(gga_x_b86_mgc_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static const FLOAT beta  = 0.00375;
  static const FLOAT gamma = 0.007;
  
  FLOAT dd, f1, f2, df1, df2, d2f1, d2f2, d3f1, d3f2;

  dd    = 1.0 + gamma*x*x;

  f1    = beta/X_FACTOR_C*x*x;
  f2    = POW(dd, 4.0/5.0);
  *f    = 1.0 + f1/f2;

  if(order < 1) return;

  df1 = beta/X_FACTOR_C*2.0*x;
  df2 = 8.0*gamma*x*f2/(5.0*dd);

  *dfdx = DFRACTION(f1, df1, f2, df2);

  if(order < 2) return;

  d2f1 = beta/X_FACTOR_C*2.0;
  d2f2 = 8.0*gamma*(5.0 + 3.0*gamma*x*x)*f2/(25.0*dd*dd);

  *d2fdx2 = D2FRACTION(f1, df1, d2f1, f2, df2, d2f2);

  if(order < 3) return;

  d3f1 = 0.0;
  d3f2 = -48.0*gamma*gamma*x*(5.0 + gamma*x*x)*f2/(125.0*dd*dd*dd);

  *d3fdx3 = D3FRACTION(f1, df1, d2f1, d3f1, f2, df2, d2f2, d3f2);
}

#define func XC(gga_x_b86_mgc_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_b86_mgc) = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  "AD Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "AD Becke, J. Chem. Phys 85, 7184 (1986)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};
