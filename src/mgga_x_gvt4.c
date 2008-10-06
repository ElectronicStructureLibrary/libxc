/*
 Copyright (C) 2008 Lara Ferrigni, Georg Madsen, M.A.L. Marques

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
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_MGGA_X_GVT4          204 /* GVT4 from Van Voorhis and Scuseria exchange */

/* calculate h and h derivatives with respect to rho, grho and tau: Equation (5) */
void XC(mgga_x_gvt4_func)(int order, FLOAT x, FLOAT z, FLOAT alpha, const FLOAT *d, 
			  FLOAT *h, FLOAT *dhdx, FLOAT *dhdz)
{
  FLOAT gam, gam2, x2, dhdgam;
  FLOAT n1, n2, n3;
  
  x2   = x*x;
  gam  = 1.0 + alpha*(x2 + z);
  gam2 = gam*gam;

  n1 = d[0];
  n2 = d[1]*x2 + d[2]*z;
  n3 = d[3]*x2*x2 + d[4]*x2*z + d[5]*z*z;

  *h = n1/gam + n2/gam2 + n3/(gam*gam2);

  if(order < 1) return;
  
  dhdgam = -n1/gam2 - 2.0*n2/(gam*gam2) - 3.0*n3/(gam2*gam2);

  *dhdx = 2.0*d[1]*x/gam2 + (4.0*d[3]*x*x2 + 2.0*d[4]*x*z)/(gam*gam2) +
    dhdgam*(2.0*alpha*x);
  *dhdz = d[2]/gam2 + (d[4]*x2 + 2.0*d[5]*z)/(gam*gam2) +
    dhdgam*alpha;
}

static void 
func(const XC(mgga_type) *pt, FLOAT x, FLOAT t, int order,
     FLOAT *f, FLOAT *dfdx, FLOAT *dfdt,
     FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  const FLOAT d[6] = {-0.9800, -0.003557, 0.006250, -0.00002354, -0.0001283, 0.0003575};
  const FLOAT alpha = 0.001867;
  const FLOAT CFermi = (3.0/5.0) * POW(6.0*M_PI*M_PI, 2.0/3.0);

  /* Eq. (14) */
  XC(mgga_x_gvt4_func)(order, x, t - CFermi, alpha, d, f, dfdx, dfdt);
 
  *f /= -X_FACTOR_C;

  if(order < 1) return;

  *dfdx /= -X_FACTOR_C;
  *dfdt /= -X_FACTOR_C;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_gvt4) = {
  XC_MGGA_X_GVT4,
  XC_EXCHANGE,
  "GVT4 (X part of VSXC)",
  XC_FAMILY_MGGA,
  "T Van Voorhis and GE Scuseria, JCP 109, 400 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
