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
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_B88          106 /* Becke 88 */

typedef struct{
  double beta;
} gga_x_b88_params;


void gga_x_b88_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;
  gga_x_b88_params *params;

  assert(p->params == NULL);

  p->params = malloc(sizeof(gga_x_b88_params));
  params = (gga_x_b88_params *) (p->params);

  /* value of beta in standard Becke 88 functional */
  gga_x_b88_set_params(p, 0.0042);
}


void gga_x_b88_end(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  assert(p->params != NULL);
  free(p->params);
}


void gga_x_b88_set_params(xc_gga_type *p, double beta)
{
  gga_x_b88_params *params;

  assert(p->params != NULL);
  params = (gga_x_b88_params *) (p->params);

  params->beta = beta;
}


static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  double f1;
  double beta;

  assert(p->params != NULL);
  beta = ((gga_x_b88_params *) (p->params))->beta;

  f1 = (1.0 + 6.0*beta*x*asinh(x));
  *f = 1.0 + beta/X_FACTOR_C*x*x/f1;
 
  *dfdx = beta/X_FACTOR_C*x*(2.0 + 6.0*beta*(x*asinh(x) - x*x/sqrt(1.0+x*x)))/(f1*f1);
  *ldfdx= beta/X_FACTOR_C;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  "A.D. Becke, Phys. Rev. A 38, 3098-3100 (1988)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_x_b88_init, 
  gga_x_b88_end, 
  NULL,
  work_gga_x
};
