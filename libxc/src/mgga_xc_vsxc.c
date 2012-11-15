/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_X_GVT4          204 /* GVT4 from Van Voorhis and Scuseria (exchange part)    */
#define XC_MGGA_C_VSXC          232 /* VSxc from Van Voorhis and Scuseria (correlation part) */

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
func(const XC(func_type) *pt, XC(work_mgga_x_params) *r)
{
  static const FLOAT abcd[6] = 
    {-9.800683e-01, -3.556788e-03, 6.250326e-03, -2.354518e-05, -1.282732e-04, 3.574822e-04};
  static const FLOAT alpha = 0.00186726;

  XC(mgga_x_gvt4_func)(r->order, r->x, 2.0*(r->t - K_FACTOR_C), alpha, abcd, &r->f, &r->dfdx, &r->dfdt);
 
  r->f /= -X_FACTOR_C;

  if(r->order < 1) return;

  r->dfdx /= -X_FACTOR_C;
  r->dfdt /= -X_FACTOR_C/2.0;
}

static void 
func_c_parallel(const XC(func_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
		FLOAT *f, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
		FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  static const FLOAT abcd[6] = {0.3271, -0.03229, -0.02942, 0.002134, -0.005452, 0.01578};
  static const FLOAT alpha = 0.005151;

  FLOAT dd, f1;
  const FLOAT tmin = 0.5e-10;

  XC(mgga_x_gvt4_func)(order, x, 2.0*(t - K_FACTOR_C), alpha, abcd, f, dfdx, dfdt);

  f1 = *f;

  dd = (t > tmin) ? 1.0 - x*x/(8.0*t) : 0.0;
  *f = dd*f1;                     /* multiply by D_sigma */

  if(order < 1) return;

  if(t > tmin){
    *dfdx  =  -x*f1/(4.0*t)   +     dd*(*dfdx);
    *dfdt  = x*x*f1/(8.0*t*t) + 2.0*dd*(*dfdt);
  }else{
    *dfdx = *dfdt = 0.0;
  }
}


static void 
func_c_opposite(const XC(func_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
		FLOAT *f, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
		FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  static const FLOAT abcd[6] = {0.7035, 0.007695, 0.05153, 0.00003394, -0.001269, 0.001296};
  static const FLOAT alpha = 0.003050;

  XC(mgga_x_gvt4_func)(order, x, 2.0*(t - 2.0*K_FACTOR_C), alpha, abcd, f, dfdx, dfdt);

  if(order < 1) return;

  *dfdt *= 2.0;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_gvt4) = {
  XC_MGGA_X_GVT4,
  XC_EXCHANGE,
  "GVT4 (X part of VSXC)",
  XC_FAMILY_MGGA,
  "T Van Voorhis and GE Scuseria, JCP 109, 400 (1998)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};

#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_vsxc) = {
  XC_MGGA_C_VSXC,
  XC_CORRELATION,
  "VSXC (correlation part)",
  XC_FAMILY_MGGA,
  "T Van Voorhis and GE Scuseria, JCP 109, 400 (1998)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU/2.0, MIN_ZETA,
  work_mgga_c_init,
  NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_c,
};
