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

#define XC_GGA_X_OPTX         110 /* Handy & Cohen OPTX 01                          */

static inline void
func(xc_gga_type *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx)
{
  static const FLOAT a = 1.05151, b = 1.43169/X_FACTOR_C, gamma = 0.006;

  FLOAT f1, u;

  f1 = 1.0 + gamma*x*x;
  u  = gamma*x*x/f1;

  *f     = a + b*u*u;
  *dfdx  = 2.0*b*u * 2.0*gamma*x/(f1*f1);
  *ldfdx = 0.0;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  "NC Handy and AJ Cohen, Mol. Phys. 99, 403 (2001)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
