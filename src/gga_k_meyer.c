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

#define XC_GGA_K_MEYER        57 /* Meyer,  Wang, and Young */

static inline void 
func(const XC(func_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT y, ll, lambda, dy, dlambda, d2lambda;

  y = X2S*x/6.0;

  ll     = LOG((1.0 + y)/ABS(1.0 - y));
  lambda = 0.5*(1.0 + (1.0 - y*y)*ll/(2.0*y));

  *f = 1.0 + lambda*x*x/(8.0*K_FACTOR_C);

  if(order < 1) return;

  dy = X2S/6.0;

  dlambda = dy*(2.0*y - (1.0 + y*y)*ll)/(4.0*y*y);

  *dfdx = (dlambda*x + 2.0*lambda)*x/(8.0*K_FACTOR_C);

  if(order < 2) return;

  d2lambda = dy*dy*(2.0*y/(y*y - 1.0) + ll)/(2.0*y*y*y);

  *d2fdx2  = (d2lambda*x*x + 4.0*dlambda*x + 2.0*lambda)/(8.0*K_FACTOR_C);
}

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_meyer) = {
  XC_GGA_K_MEYER,
  XC_KINETIC,
  "Meyer,  Wang, and Young",
  XC_FAMILY_GGA,
  "A Meyer, GC Wang and WH Young, Z. Naturforsch. A 31, 898-903 (1976)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL,
  NULL, NULL,
  work_gga_k,
  NULL
};
