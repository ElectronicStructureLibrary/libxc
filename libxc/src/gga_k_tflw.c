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

#define XC_GGA_K_TFLW          201 /* Thomas-Fermi-von Weizsaecker-like functionals */

typedef struct{
  FLOAT lambda;
} gga_k_tflw_params;


static void 
gga_k_tflw_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_k_tflw_params));

  /* possible values of lambda can be found in the reference. 
     However, I do not have access to the article ;(  */

  XC(gga_k_tflw_set_params_)(p, 1.0);
}

void 
XC(gga_k_tflw_set_params)(XC(func_type) *p, FLOAT lambda)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_k_tflw_set_params_)(p->gga, lambda);
}


void 
XC(gga_k_tflw_set_params_)(XC(gga_type) *p, FLOAT lambda)
{
  gga_k_tflw_params *params;

  assert(p->params != NULL);
  params = (gga_k_tflw_params *) (p->params);

  params->lambda  = lambda;
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  FLOAT lambda, ss;

  assert(p->params != NULL);
  lambda  = ((gga_k_tflw_params *) (p->params))->lambda;
  lambda *= 5.0/3.0;

  ss  = X2S*x;

  *f = 1.0 + lambda*ss*ss;

  if(order < 1) return;

  *dfdx = X2S*2.0*lambda*ss;
  *ldfdx= X2S*X2S*lambda;
  
  if(order < 2) return;

  *d2fdx2 = X2S*X2S*2.0*lambda;
}

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_tflw) = {
  XC_GGA_K_TFLW,
  XC_KINETIC,
  "Thomas-Fermi-von Weizsaecker-like functional (TFlambdaW)",
  XC_FAMILY_GGA,
  "E. V. Ludeña and V. V. Karasiev, in 'Reviews of Modern Quantum Chemistry: A Celebration of the Contributions of Robert Parr', edited by K. D. Sen (World Scientific, Singapore, 2002), pp. 612-665.",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_x
};
