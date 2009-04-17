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

#define XC_GGA_X_RPBE  117 /* Hammer, Hansen & Norskov (PBE-like) */

/* RPBE: see PBE for more details */
static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT kappa = 0.8040;
  static const FLOAT mu = 0.00361218645365094697;

  FLOAT f0, df0, d2f0;

  f0 = exp(-mu*x*x/kappa);
  *f = 1.0 + kappa*(1.0 - f0);

  if(order < 1) return;

  df0 = -2.0*x*mu/kappa*f0;
  
  *dfdx  = -kappa*df0;
  *ldfdx = mu;

  if(order < 2) return;

  d2f0    = -2.0*mu/kappa*f0*(1.0 - 2.0*x*x*mu/kappa);
  *d2fdx2 = -kappa*d2f0;
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_rpbe) = {
  XC_GGA_X_RPBE,
  XC_EXCHANGE,
  "Hammer, Hansen, and Nørskov",
  XC_FAMILY_GGA,
  "B Hammer, LB Hansen and JK Nørskov, Phys. Rev. B 59, 7413 (1999)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};
