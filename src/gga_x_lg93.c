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
func(xc_gga_type *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx)
{
  static const FLOAT x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  static const FLOAT ad = 1e-8, a4 = 29.790, a6 = 22.417;
  static const FLOAT a8 = 12.119, a10 = 1570.1, a12 = 55.944;
  static const FLOAT a2 = 4.94113918475214219939; /* (ad + 0.1234)/b, b = 0.024974 */

  FLOAT ss, ss2, ss4, ss6, ss8, ss10;
  FLOAT f1, f2, f3;

  ss  = x2s*x;    ss2  = ss*ss;
  ss4 = ss2*ss2;  ss6  = ss4*ss2;
  ss8 = ss6*ss2;  ss10 = ss8*ss2;

  f1 = 1.0 + a2*ss2 + a4*ss4 + a6*ss6 + a8*ss8 + a10*ss10 + a12*ss2*ss10;
  f2 = 1.0 + ad*ss2;

  *f = f1/f2;

  f3 = 2.0*ss*(a2 + 2.0*a4*ss2 + 3.0*a6*ss4 + 4.0*a8*ss6 + 5.0*a10*ss8 + 6.0*a12*ss10);
  *dfdx  = x2s*(f3*f2 - 2.0*ss*ad*f1)/(f2*f2);
  *ldfdx = x2s*x2s*(a2 - ad);
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_lg93 = {
  XC_GGA_X_LG93,
  XC_EXCHANGE,
  "Lacks & Gordon 93",
  XC_FAMILY_GGA,
  "DJ Lacks and RG Gordon, Phys. Rev. A 47, 4681 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

