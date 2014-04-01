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

#define XC_GGA_X_DK87_R1      111 /* dePristo & Kress 87 (version R1)               */
#define XC_GGA_X_DK87_R2      112 /* dePristo & Kress 87 (version R2)               */

void XC(gga_x_dk87_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static const FLOAT a1[2] = {0.861504, 0.861213}, 
    b1[2] = {0.044286, 0.042076}, alpha[2] = {1.0, 0.98};
  static const FLOAT betag = 0.00132326681668994855/X_FACTOR_C; /* 7/(432*pi*(6*pi^2)^(1/3)) */
  
  FLOAT f0, f1, f2, df1, df2, d2f1, d2f2, d3f1;
  int func;

  switch(p->info->number){
  case XC_GGA_X_DK87_R2: func = 1; break;
  default:               func = 0; /* XC_GGA_X_DK87_R1 */
  }

  f0 = a1[func]*POW(x, alpha[func]);
  f1 = betag*x*x*(1.0 + f0);
  f2 = 1.0 + b1[func]*x*x;
  
  *f     = 1.0 + f1/f2;

  if(order < 1) return;

  df1 = betag*x*(2.0 + (2.0 + alpha[func])*f0);
  df2 = 2.0*b1[func]*x;

  *dfdx  = DFRACTION(f1, df1, f2, df2);
  
  if(order < 2) return;

  d2f1 = betag*(2.0 + (1.0 + alpha[func])*(2.0 + alpha[func])*f0);
  d2f2 = 2.0*b1[func];

  *d2fdx2 = D2FRACTION(f1, df1, d2f1, f2, df2, d2f2);

  if(order < 3) return;

  d3f1 = betag*alpha[func]*(1.0 + alpha[func])*(2.0 + alpha[func])*f0/x;
  
  *d3fdx3 = D3FRACTION(f1, df1, d2f1, d3f1, f2, df2, d2f2, 0.0);
}

#define func XC(gga_x_dk87_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_dk87_r1) = {
  XC_GGA_X_DK87_R1,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R1",
  XC_FAMILY_GGA,
  "AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24, 1e-24, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_dk87_r2) = {
  XC_GGA_X_DK87_R2,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R2",
  XC_FAMILY_GGA,
  "AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};
