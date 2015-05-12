/*
 Copyright (C) 2015 Susi Lehtola

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

#define XC_MGGA_X_MVS          257 /* MVS exchange of Sun, Perdew, and Ruzsinszky */

static void
func_fa(int order, FLOAT a, 
	FLOAT *f, FLOAT *dfda, FLOAT *d2fda2)
{
  const FLOAT e1=-1.6665;
  const FLOAT c1= 0.7438;
  
  FLOAT a2, a4, aux, num, den, bareden, dnum, dden, d2den;
  FLOAT e12;

  e12 = e1*e1;
  
  a2 = a*a;
  a4 = a2*a2;
  
  num = 1 - a;

  aux = 1 + e1*a2;
  /* Use double square root instead of pow, this is faster */
  bareden = aux*aux + c1*a4;
  den=SQRT(SQRT(bareden));
  
  *f = num/den;

  if(order < 1) return;              /* And now the derivatives */

  dnum  = -1.0;
  dden  = a*(e1 + a2*(e12 + c1))*(den/bareden);

  *dfda = (dnum*den - num*dden)/(den*den);

  if(order < 2) return;

  d2den   = (e1*(e12+c1)*a4 + (2.0*e12 + 3.0*c1)*a2 + e1)*den/(bareden*bareden);

  *d2fda2 = (-num*d2den*den - 2.0*dden*(dnum*den - dden*num))/(den*den*den);
}

static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  FLOAT ss;
  FLOAT a, dadt, dadx, d2adx2, fa, dfada, d2fada2;
  FLOAT p, dpdx, d2pdx2;
  FLOAT opbp, opbprt, dfda, dfdp, d2fda2, d2fdp2, d2fdadp;

  const FLOAT k0 = 0.174;
  const FLOAT b = 0.0233;

  ss = X2S*r->x;
  p = ss*ss;

  a = (r->t - r->x*r->x/8.0)/K_FACTOR_C;

  /* Calculate f(alpha) */
  func_fa(r->order, a, &fa, &dfada, &d2fada2);

  /* Calculate helpers */
  opbp = 1.0 + b*p*p;
  opbprt = POW(opbp,1.0/8.0);

  /* Functional value is */
  r->f = (1.0 + k0*fa)/opbprt;

  if(r->order < 1) return;

  dpdx = 2.0*ss*X2S;
  dadx = -2.0*r->x/(8.0*K_FACTOR_C);
  dadt = 1.0/K_FACTOR_C;

  dfda = k0/opbprt*dfada;
  dfdp = -1.0/8.0*(1 + k0*fa)/opbp*opbprt;

  r->dfdx = dfda*dadx + dfdp*dpdx;
  r->dfdt = dfda*dadt;
  r->dfdu = 0.0;

  if(r->order < 2) return;

  d2pdx2 = 2.0*X2S*X2S;
  d2adx2 = -2.0/(8.0*K_FACTOR_C);

  d2fda2 =  k0/opbprt*d2fada2;
  d2fdp2 = 7.0/64.0*(1 + k0*fa)/(opbp*opbp)*opbprt;
  d2fdadp = -1.0/8.0*k0*dfada/opbp*opbprt;
  
  r->d2fdx2 = dfdp*d2pdx2 + d2fdp2*dpdx*dpdx + dfda*d2adx2 + d2fda2*dadx*dadx;
  r->d2fdt2 = d2fda2*dadt*dadt;
  r->d2fdxt = d2fdadp*dadt*dadx;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_mvs) = {
  XC_MGGA_X_MVS,
  XC_EXCHANGE,
  "MVS exchange of Sun, Perdew, and Ruzsinszky",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_685, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_DEVELOPMENT,
  1e-32, 1e-32, 1e-32, 1e-32,
  NULL, NULL, NULL, NULL,
  work_mgga_x,
};
