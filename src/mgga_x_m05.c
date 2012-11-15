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

#define XC_MGGA_X_M05          214 /* M05 functional of Minnesota */
#define XC_MGGA_X_M05_2X       215 /* M05-2X functional of Minnesota */

static const FLOAT a_m05[12] = 
  {1.0, 0.08151, -0.43956, -3.22422, 2.01819, 8.79431, -0.00295,
   9.82029, -4.82351, -48.17574, 3.64802, 34.02248};

static const FLOAT a_m05_2x[12] =
  {1.0, -0.56833, -1.30057, 5.50070, 9.06402, -32.21075, -23.73298,
   70.22996, 29.88614, -60.25778, -13.22205, 15.23694};

typedef struct{
  int n;
  const FLOAT *a;
} mgga_x_m05_params;


static void
mgga_x_m05_init(XC(func_type) *p)
{
  mgga_x_m05_params *params;

  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_X_PBE, p->nspin);

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m05_params));
  params = (mgga_x_m05_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_M05: 
    params->n = 12;
    params->a = a_m05;
    break;
  case XC_MGGA_X_M05_2X:
    params->n = 12;
    params->a = a_m05_2x;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_tpss\n");
    exit(1);
  }
}


static void 
func(const XC(func_type) *pt, XC(work_mgga_x_params) *r)
{
  mgga_x_m05_params *params;

  FLOAT e_f, e_dfdx, e_d2fdx2;
  FLOAT w, w_den, wi, factor;
  FLOAT dwdt, dfactordw;
  int i;

  assert(pt != NULL && pt->params != NULL);
  params = (mgga_x_m05_params *) (pt->params);
  
  XC(gga_x_pbe_enhance)(pt->func_aux[0], r->order, r->x, &e_f, &e_dfdx, &e_d2fdx2);

  w_den = K_FACTOR_C + r->t;
  w = (K_FACTOR_C - r->t)/w_den;

  factor = 0.0;
  wi     = 1.0;

  for(i=0; i<params->n; i++){
    factor += params->a[i]*wi;
    wi *= w;
  }

  r->f = e_f*factor;

  if(r->order < 1) return;

  dwdt = -2.0*K_FACTOR_C/(w_den*w_den);

  dfactordw = 0.0;
  wi        = 1.0;
  for(i=1; i<params->n; i++){
    dfactordw += i*params->a[i]*wi;
    wi *= w;
  }

  r->dfdx = e_dfdx*factor;
  r->dfdt = e_f*dfactordw*dwdt;
  r->dfdu = 0.0;

  if(r->order < 2) return;

}


#include "work_mgga_x.c"


XC(func_info_type) XC(func_info_mgga_x_m05) = {
  XC_MGGA_X_M05,
  XC_EXCHANGE,
  "M05 functional of Minnesota",
  XC_FAMILY_MGGA,
  "Y Zhao, NE Schultz, and DG Truhlar, J. Chem. Phys. 123, 161103 (2005)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m05_init,
  NULL, NULL, NULL,
  work_mgga_x,
};


XC(func_info_type) XC(func_info_mgga_x_m05_2x) = {
  XC_MGGA_X_M05_2X,
  XC_EXCHANGE,
  "M05-2X functional of Minnesota",
  XC_FAMILY_MGGA,
  "Y Zhao, NE Schultz, and DG Truhlar, J. Chem. Theory Comput. 2, 364 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m05_init,
  NULL, NULL, NULL,
  work_mgga_x,
};
