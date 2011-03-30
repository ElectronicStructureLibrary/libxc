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

#define XC_GGA_K_DK          516 /* DePristo and Kress */

static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT aa[5] = {1.0,  0.95, 14.28111, -19.57962, 26.64711};
  static const FLOAT bb[4] = {1.0, -0.05,  9.99802,   2.96085};

  FLOAT ss, ss2, num, denom, dnum, ddenom, d2num, d2denom;

  ss  = x/(72.0*X_FACTOR_C);
  ss2 = ss*ss;

  num   = aa[0] + aa[1]*ss + aa[2]*ss2 + aa[3]*ss*ss2 + 9.0*bb[3]*ss2*ss2;
  denom = bb[0] + bb[1]*ss + bb[2]*ss2 + bb[3]*ss*ss2;

  *f = num/denom;

  if(order < 1) return;

  dnum   = aa[1] + 2.0*aa[2]*ss + 3.0*aa[3]*ss2 + 4.0*9.0*bb[3]*ss*ss2;
  ddenom = bb[1] + 2.0*bb[2]*ss + 3.0*bb[3]*ss2;

  *dfdx  = (dnum*denom - num*ddenom)/(denom*denom);
  *dfdx /= 72.0*X_FACTOR_C;
  
  if(order < 2) return;

  d2num   = 2.0*aa[2] + 3.0*2.0*aa[3]*ss + 4.0*3.0*9.0*bb[3]*ss2;
  d2denom = 2.0*bb[2] + 3.0*2.0*bb[3]*ss;

  *d2fdx2  = ((d2num*denom - num*d2denom)*denom - 2.0*ddenom*(dnum*denom - ddenom*num))/(denom*denom*denom);
  *d2fdx2 /= 72.0*X_FACTOR_C * 72.0*X_FACTOR_C;
}

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_dk) = {
  XC_GGA_K_DK,
  XC_KINETIC,
  "DePristo and Kress",
  XC_FAMILY_GGA,
  "AE DePristo and JD Kress, Phys. Rev. A 35, 438-441 (1987)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  NULL, NULL, NULL,
  work_gga_k
};
