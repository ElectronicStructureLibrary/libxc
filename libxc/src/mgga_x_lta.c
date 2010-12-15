/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

#include <stdlib.h>
#include <assert.h>
#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_LTA          201 /* Local tau approximation of Ernzerhof & Scuseria */

static void 
func(const XC(mgga_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
     FLOAT *f, FLOAT *vrho0, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
     FLOAT *d2fdx2, FLOAT *d2fdt2, FLOAT *d2fdu2, FLOAT *d2fdxt, FLOAT *d2fdxu, FLOAT *d2fdtu)
{
  /* POW(10.0/(3.0*POW(6.0*M_PI*M_PI, 2.0/3.0)), 4.0/5.0) = (2/C_F)^(4/5) */
  const FLOAT a1 = 0.297163728291293581339216378935;

  t  = t/2.0; /* we use a different definition of t */
  *f = a1*POW(t, 4.0/5.0);

  if(order < 1) return;
  
  *dfdx = 0.0;
  *dfdt = (t > 1e-10) ? a1*4.0/5.0*POW(t, -1.0/5.0)/2.0 : 0.0;
  *dfdu = 0.0;

  if(order < 2) return;
  
  *d2fdx2 = 0.0;
  *d2fdt2 = (t > 1e-10) ? -a1*4.0/25.0*POW(t, -6.0/5.0)/4.0 : 0.0;
  *d2fdu2 = 0.0;
  *d2fdxt = 0.0;
  *d2fdxu = 0.0;
  *d2fdtu = 0.0;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_lta) = {
  XC_MGGA_X_LTA,
  XC_EXCHANGE,
  "Local tau approximation",
  XC_FAMILY_MGGA,
  "M Ernzerhof and G Scuseria, J. Chem. Phys. 111, 911 (1999)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
