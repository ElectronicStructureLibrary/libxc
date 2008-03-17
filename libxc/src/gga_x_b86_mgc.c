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

#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalfa,beta,gamma (with mod. grad. correction) */

static inline void
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT beta  = 0.00375;
  static const FLOAT gamma = 0.007;
  
  FLOAT f1;

  f1    = (1.0 + gamma*x*x);
  *f    = 1.0 + beta/X_FACTOR_C*x*x/POW(f1, 4.0/5.0);

  *dfdx = beta/X_FACTOR_C*2.0*x*(5.0 + gamma*x*x)/(5.0*POW(f1, 9.0/5.0));
  *ldfdx= beta/X_FACTOR_C;
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_b86_mgc) = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  "AD Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "AD Becke, J. Chem. Phys 85, 7184 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
