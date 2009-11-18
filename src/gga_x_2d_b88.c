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

#define XC_GGA_X_2D_B88        127 /* Becke 88 in 2D */

typedef struct{
  FLOAT beta;
} gga_x_2d_b88_params;

void XC(gga_x_2d_b88_set_params_)(XC(gga_type) *p, FLOAT beta);

static void gga_x_2d_b88_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  gga_x_2d_b88_params *params;

  assert(p->params == NULL);

  p->params = malloc(sizeof(gga_x_2d_b88_params));
  params = (gga_x_2d_b88_params *) (p->params);

  /* value of beta in standard Becke 88 2D functional */
  XC(gga_x_2d_b88_set_params_)(p, 0.018641);
}


static void gga_x_2d_b88_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}


void XC(gga_x_2d_b88_set_params)(XC(func_type) *p, FLOAT beta)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_x_2d_b88_set_params_)(p->gga, beta);
}


void XC(gga_x_2d_b88_set_params_)(XC(gga_type) *p, FLOAT beta)
{
  gga_x_2d_b88_params *params;

  assert(p->params != NULL);
  params = (gga_x_2d_b88_params *) (p->params);

  params->beta = beta;
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  FLOAT f1, f2, df1, df2, d2f1, d2f2;
  FLOAT beta, csi;

  assert(p->params != NULL);
  beta = ((gga_x_2d_b88_params *) (p->params))->beta;
  csi  = 8.0; /* for harmonic potentials */

  f1 = beta/X_FACTOR_2D_C*x*x;
  f2 = 1.0 + csi*beta*x*asinh(x);
  *f = 1.0 + f1/f2;

  if(order < 1) return;

  df1 = 2.0*beta/X_FACTOR_2D_C*x;
  df2 = csi*beta*(asinh(x) + x/sqrt(1.0 + x*x));

  *dfdx = (df1*f2 - f1*df2)/(f2*f2);
  *ldfdx= beta/X_FACTOR_2D_C;

  if(order < 2) return;

  d2f1 = 2.0*beta/X_FACTOR_2D_C;
  d2f2 = csi*beta*(2.0 + x*x)/pow(1.0 + x*x, 3.0/2.0);

  *d2fdx2 = (2.0*f1*df2*df2 + d2f1*f2*f2 - f2*(2.0*df1*df2 + f1*d2f2))/(f2*f2*f2);
}

#define XC_DIMENSIONS 2
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_2d_b88) = {
  XC_GGA_X_2D_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  "G Vilhena, MAL Marques, unpublished\n"
  "AD Becke, Phys. Rev. A 38, 3098 (1988)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  gga_x_2d_b88_init, 
  gga_x_2d_b88_end, 
  NULL,
  work_gga_x
};
