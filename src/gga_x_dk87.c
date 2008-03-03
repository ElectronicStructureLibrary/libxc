/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_DK87_R1      111 /* dePristo & Kress 87 (version R1)               */
#define XC_GGA_X_DK87_R2      112 /* dePristo & Kress 87 (version R2)               */

static inline void 
func(xc_gga_type *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx)
{
  static const FLOAT a1[2] = {0.861504, 0.861213}, 
    b1[2] = {0.044286, 0.042076}, alpha[2] = {1.0, 0.98};
  static const FLOAT betag = 0.00132326681668994855/X_FACTOR_C; /* 7/(432*pi*(6*pi^2)^(1/3)) */
  
  FLOAT f0, f1, f2;
  int func;

  switch(p->info->number){
  case XC_GGA_X_DK87_R2: func = 1; break;
  default:               func = 0; /* XC_GGA_X_DK87_R1 */
  }

  f0 = a1[func]*POW(x, alpha[func]);
  f1 = 1.0 + f0;
  f2 = 1.0 + b1[func]*x*x;
  
  *f     = 1.0 + betag*x*x*f1/f2;
  *dfdx  = betag*(2.0*x*f1/f2 + x*(alpha[func]*f0*f2 - 2.0*b1[func]*x*x*f1)/(f2*f2));
  *ldfdx = betag;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_dk87_r1 = {
  XC_GGA_X_DK87_R1,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R1",
  XC_FAMILY_GGA,
  "AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

const xc_func_info_type func_info_gga_x_dk87_r2 = {
  XC_GGA_X_DK87_R2,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R2",
  XC_FAMILY_GGA,
  "AE DePristo and JD Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
