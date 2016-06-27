/*
 Copyright (C) 2014 Jess Wellendorff, M.A.L. Marques

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

#define XC_GGA_X_BEEFVDW          253 /* BEEF-vdW exchange */
#define XC_GGA_XC_BEEFVDW         256 /* BEEF-vdW exchange-correlation */

static const FLOAT coefs_beefvdw[30] = {
   1.516501714e0,   4.413532099e-1, -9.182135241e-2, -2.352754331e-2,
   3.418828455e-2,  2.411870076e-3, -1.416381352e-2,  6.975895581e-4,
   9.859205137e-3, -6.737855051e-3, -1.573330824e-3,  5.036146253e-3,
  -2.569472453e-3, -9.874953976e-4,  2.033722895e-3, -8.018718848e-4,
  -6.688078723e-4,  1.030936331e-3, -3.673838660e-4, -4.213635394e-4,
   5.761607992e-4, -8.346503735e-5, -4.458447585e-4,  4.601290092e-4,
  -5.231775398e-6, -4.239570471e-4,  3.750190679e-4,  2.114938125e-5,
  -1.904911565e-4,  7.384362421e-5
};

typedef struct{
  int legorder;
  const FLOAT *coefs;
} gga_x_beef_params;

static void
gga_x_beef_init(XC(func_type) *p)
{
  gga_x_beef_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_beef_params));
  params = (gga_x_beef_params *)p->params;

  switch(p->info->number){
  case XC_GGA_X_BEEFVDW:
    params->legorder = 30;
    params->coefs = coefs_beefvdw;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_beef\n");
    exit(1);
  }
}

static inline void
transformations(int order,FLOAT p, FLOAT *xi, FLOAT *dxidp)
{
  FLOAT tmp0, num, den, dnum, dden;
  const FLOAT k = 4.0;

  /* reduced density gradient in transformation t1(s) */
  tmp0   = k + p;
  *xi    = 2.0 * p / tmp0 - 1.0;

  if(order < 1) return;

  *dxidp = 2.0 * k / (tmp0*tmp0);
}

static inline void
legf(int order, int legorder, const FLOAT *coefs, FLOAT xi, 
     FLOAT *f, FLOAT *dfdxi)
{
  int i, j;
  FLOAT Li[legorder], dLi[legorder];

  /* product exchange enhancement factor and derivatives */

  /* initializing */
  Li[0] = 1.0; Li[1] = xi;

  /* recursively building polynomials and their derivatives */
  for(i = 2; i < legorder; i++){
    Li[i] = 2.0 * xi * Li[i-1] - Li[i-2] - (xi * Li[i-1] - Li[i-2])/i;
  }

  /* building enhancement factor and derivatives */
  *f = 0.0;
  for(j = 0; j < legorder; j++){
    *f += coefs[j]*Li[j];
  }

  if(order < 1) return;

  dLi[0] = 0.0; dLi[1] = 1.0;
  
  for(i = 2; i < legorder; i++){
    dLi[i] = i * Li[i-1] + xi * dLi[i-1];
  }
  
  *dfdxi = 0.0;
  for(j = 0; j < legorder; j++){
      *dfdxi += coefs[j]*dLi[j];
  }  
}

static void 
func(const XC(func_type) *pt, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT p, dpdx;
  FLOAT xi, dxidp, fx, dfdxi;
  gga_x_beef_params *params;

  params = (gga_x_beef_params *)pt->params;
  assert(params != NULL);

  p = X2S*X2S*x*x;

  transformations(order, p, &xi, &dxidp);
  legf(order, params->legorder, params->coefs, xi, &fx, &dfdxi);

  *f = fx;

  if(order < 1) return;

  dpdx = 2.0*X2S*X2S*x;

  *dfdx = dfdxi*dxidp*dpdx;
}


#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_beefvdw) = {
  XC_GGA_X_BEEFVDW,
  XC_EXCHANGE,
  "BEEF-vdW exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_x_beef_init,
  NULL, NULL, work_gga_x, NULL,
};


void
gga_xc_beefvdw_init(XC(func_type) *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_BEEFVDW, XC_LDA_C_PW_MOD, XC_GGA_C_PBE};
  static FLOAT funcs_coef[3] = {1.0, 0.6001664769, 1.0 - 0.6001664769};

  XC(mix_init)(p, 3, funcs_id, funcs_coef);
}

const XC(func_info_type) XC(func_info_gga_xc_beefvdw) = {
  XC_GGA_XC_BEEFVDW,
  XC_EXCHANGE_CORRELATION,
  "BEEF-vdW exchange-correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_xc_beefvdw_init,
  NULL, NULL, NULL, NULL,
};
