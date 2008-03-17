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
#include "util.h"

#define XC_GGA_X_LG93  113 /* Lacks & Gordon 93 */

static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  static const FLOAT ad = 1e-8, a4 = 29.790, a6 = 22.417;
  static const FLOAT a8 = 12.119, a10 = 1570.1, a12 = 55.944;
  static const FLOAT a2 = 4.94113918475214219939; /* (ad + 0.1234)/b, b = 0.024974 */

  FLOAT ss, ss2, ss4, ss6, ss8, ss10;
  FLOAT f1, f2, df1, df2, d2f1, d2f2;

  ss  = x2s*x;    ss2  = ss*ss;
  ss4 = ss2*ss2;  ss6  = ss4*ss2;
  ss8 = ss6*ss2;  ss10 = ss8*ss2;

  f1 = 1.0 + a2*ss2 + a4*ss4 + a6*ss6 + a8*ss8 + a10*ss10 + a12*ss2*ss10;
  f2 = 1.0 + ad*ss2;

  *f = f1/f2;

  /* now come the first derivatives */
  if(dfdx==NULL && d2fdx2==NULL) return; /* nothing else to do */

  df1 = 2*ss*(a2 + 2*a4*ss2 + 3*a6*ss4 + 4*a8*ss6 + 5*a10*ss8 + 6*a12*ss10);
  df2 = 2*ss*ad;

  if(dfdx!=NULL){
    *dfdx  = x2s*(df1*f2 - f1*df2)/(f2*f2);
    *ldfdx = x2s*x2s*(a2 - ad);
  }

  if(d2fdx2==NULL) return; /* nothing else to do */
  
  d2f1 = 2*1*a2 + 4*3*a4*ss2 + 6*5*a6*ss4 + 8*7*a8*ss6 + 10*9*a10*ss8 + 12*11*a12*ss10;
  d2f2 = 2*ad;

  *d2fdx2 = x2s*x2s*(2*f1*df2*df2 + d2f1*f2*f2 - f2*(2*df1*df2 + f1*d2f2))/(f2*f2*f2);
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_lg93) = {
  XC_GGA_X_LG93,
  XC_EXCHANGE,
  "Lacks & Gordon 93",
  XC_FAMILY_GGA,
  "DJ Lacks and RG Gordon, Phys. Rev. A 47, 4681 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};

