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

#define XC_GGA_X_2D_B86_MGC      124 /* Becke 86 MGC for 2D systems */

static inline void
func(const XC(gga_type) *p, int order, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT beta=0.003317, gam=0.008323;

  FLOAT dd, ddp, f1, f2, df1, df2, d2f1, d2f2;

  dd    = 1.0 + gam*x*x;

  f1    = beta/X_FACTOR_C*x*x;
  f2    = POW(dd, 3.0/4.0);
  *f    = 1.0 + f1/f2;

  if(order < 1) return; /* nothing else to do */

  df1 = beta/X_FACTOR_C*2.0*x;
  ddp = gam*2.0*3.0/4.0*f2/dd;
  df2 = ddp*x;

  *dfdx  = (df1*f2 - f1*df2)/(f2*f2);
  *ldfdx = beta/X_FACTOR_C;

  if(order < 2) return; /* nothing else to do */

  d2f1 = beta/X_FACTOR_C*2.0;
  d2f2 = ddp*(1.0 - 2.0/4.0*gam*x*x/dd);

  *d2fdx2 = (2.0*f1*df2*df2 + d2f1*f2*f2 - f2*(2.0*df1*df2 + f1*d2f2))/(f2*f2*f2);
}

#define XC_DIMENSIONS 2
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_2d_b86_mgc) = {
  XC_GGA_X_2D_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction for 2D",
  XC_FAMILY_GGA,
  "S Pittalis, E Rasanen, JG Vilhena, and MAL Marques, Phys. Rev. A 79, 012503 (2009)\n"
  "AD Becke, J. Chem. Phys 85, 7184 (1986)",
  XC_FLAGS_2D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  NULL, NULL, NULL,
  work_gga_x
};
