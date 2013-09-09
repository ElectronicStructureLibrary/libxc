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

#define XC_GGA_X_LV_RPW86 58 /* Berland and Hyldgaard */

static void
gga_x_lv_rpw86_init(XC(func_type) *p)
{
  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_X_RPW86, p->nspin);
}


void XC(gga_x_lv_rpw86_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2)
{
  const FLOAT alpha=0.02178, beta=1.15, muLV=0.8491/9.0;
  FLOAT ss, ss2, ss4, ss6, a, da, d2a;
  FLOAT aux1, aux2;

  ss  = X2S*x;
  ss2 = ss*ss;
  ss4 = ss2*ss2;
  ss6 = ss2*ss4;

  aux1 = 1.0 + alpha*ss6;
  aux2 = beta + alpha*ss6;

  XC(gga_x_pw86_enhance)(p->func_aux[0], order, x, &a, &da, &d2a);

  *f = (1.0 + muLV*ss2)/aux1 + alpha*ss6*a/aux2;

  if(order < 1) return;
  da /= X2S;

  *dfdx = ss*(2.0*muLV - 6.0*alpha*ss4 - 4.0*alpha*muLV*ss6)/(aux1*aux1) +
    alpha*ss*ss4*(6.0*beta*a + ss*aux2*da)/(aux2*aux2);

  *dfdx *= X2S;

  if(order < 2) return;
  d2a /= X2S*X2S;

  *d2fdx2 = (6.0*alpha*ss4*(7.0*alpha*ss6 - 5.0) + 2.0*muLV*(1.0 + 5.0*alpha*ss6*(2.0*alpha*ss6 - 5.0)))/(aux1*aux1*aux1) +
    alpha*ss4*(6.0*beta*a*(5.0*beta - 7.0*alpha*ss6) + ss*aux2*(12.0*beta*da + ss*aux2*d2a))/(aux2*aux2*aux2);

  *d2fdx2 *= X2S*X2S;
}

#define func XC(gga_x_lv_rpw86_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_lv_rpw86) = {
  XC_GGA_X_LV_RPW86,
  XC_EXCHANGE,
  "Berland and Hyldgaard",
  XC_FAMILY_GGA,
  "K Berland and P Hyldgaard, arXiv:1309.1756 [cond-mat.mtrl-sci] (2013)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_lv_rpw86_init,
  NULL, NULL, 
  work_gga_x,
  NULL
};
