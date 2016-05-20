/*
 Copyright (C) 2006-2014 L. Talirz, M.A.L. Marques

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

#define XC_GGA_X_B86          103 /* Becke 86 Xalpha,beta,gamma                      */
#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_B86_R         41 /* Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */

typedef struct{
  FLOAT beta, gamma, omega;
} gga_x_b86_params;


static void 
gga_x_b86_init(XC(func_type) *p)
{
  FLOAT mu, kappa;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_b86_params));

  /* value of beta and gamma in Becke 86 functional */
  switch(p->info->number){
  case XC_GGA_X_B86:
    p->func = 0; 
    XC(gga_x_b86_set_params)(p, 0.0036/X_FACTOR_C, 0.004, 1.0);
    break;
  case XC_GGA_X_B86_MGC:
    p->func = 2; 
    XC(gga_x_b86_set_params)(p, 0.00375/X_FACTOR_C, 0.007, 4.0/5.0);
    break;
  case XC_GGA_X_B86_R:
    p->func = 2; 
    mu = 10.0/81.0;
    kappa = 0.7114;
    XC(gga_x_b86_set_params)(p, mu*X2S*X2S, mu*X2S*X2S/kappa, 4.0/5.0);
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_b86\n");
    exit(1);
  }
}


void 
XC(gga_x_b86_set_params)(XC(func_type) *p, FLOAT beta, FLOAT gamma, FLOAT omega)
{
  gga_x_b86_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b86_params *) (p->params);

  params->beta  = beta;
  params->gamma = gamma;
  params->omega = omega;
}


void XC(gga_x_b86_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT beta, gamma, omega;
  FLOAT dd, ddd, d2dd;
  FLOAT f1, f2, df1, df2, d2f1, d2f2, d3f2;

  assert(p->params != NULL);
  beta  = ((gga_x_b86_params *) (p->params))->beta;
  gamma = ((gga_x_b86_params *) (p->params))->gamma;
  omega = ((gga_x_b86_params *) (p->params))->omega;

  dd = 1.0 + gamma*x*x;
  f1 = beta*x*x;
  f2 = POW(dd, omega);

  *f = 1.0 + f1/f2;
  
  if(order < 1) return;

  ddd = 2.0*gamma*x;
  df1 = 2.0*beta *x;
  df2 = omega*ddd*f2/dd;

  *dfdx  = DFRACTION(f1, df1, f2, df2);

  if(order < 2) return;

  d2dd = 2.0*gamma;
  d2f1 = 2.0*beta;
  d2f2 = omega*f2/(dd*dd)*(d2dd*dd + (omega - 1.0)*ddd*ddd);

  *d2fdx2 = D2FRACTION(f1, df1, d2f1, f2, df2, d2f2);

  if(order < 3) return;

  d3f2 = omega*(omega - 1.0)*ddd*f2/(dd*dd*dd)*(3.0*d2dd*dd + (omega - 2.0)*ddd*ddd);

  *d3fdx3 = D3FRACTION(f1, df1, d2f1, 0.0, f2, df2, d2f2, d3f2);
}

#define func XC(gga_x_b86_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_b86) = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-23, 0.0, 1e-32,
  gga_x_b86_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_b86_mgc) = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_b86_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_b86_r) = {
  XC_GGA_X_B86_R,
  XC_EXCHANGE,
  "Revised Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Hamada2014_121103, &xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_b86_init,
  NULL, NULL,
  work_gga_x,
  NULL
};
