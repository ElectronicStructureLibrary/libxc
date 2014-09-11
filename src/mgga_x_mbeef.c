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

#define XC_MGGA_X_MBEEF          249 /* mBEEF exchange */
#define XC_MGGA_X_MBEEFVDW       250 /* mBEEF-vdW exchange */

static const FLOAT coefs_mbeef[64] = {
         1.18029330e+00,   8.53027860e-03,  -1.02312143e-01,
         6.85757490e-02,  -6.61294786e-03,  -2.84176163e-02,
         5.54283363e-03,   3.95434277e-03,  -1.98479086e-03,
         1.00339208e-01,  -4.34643460e-02,  -1.82177954e-02,
         1.62638575e-02,  -8.84148272e-03,  -9.57417512e-03,
         9.40675747e-03,   6.37590839e-03,  -8.79090772e-03,
        -1.50103636e-02,   2.80678872e-02,  -1.82911291e-02,
        -1.88495102e-02,   1.69805915e-07,  -2.76524680e-07,
         1.44642135e-03,  -3.03347141e-03,   2.93253041e-03,
        -8.45508103e-03,   6.31891628e-03,  -8.96771404e-03,
        -2.65114646e-08,   5.05920757e-08,   6.65511484e-04,
         1.19130546e-03,   1.82906057e-03,   3.39308972e-03,
        -7.90811707e-08,   1.62238741e-07,  -4.16393106e-08,
         5.54588743e-08,  -1.16063796e-04,   8.22139896e-04,
        -3.51041030e-04,   8.96739466e-04,   2.09603871e-08,
        -3.76702959e-08,   2.36391411e-08,  -3.38128188e-08,
        -5.54173599e-06,  -5.14204676e-05,   6.68980219e-09,
        -2.16860568e-08,   9.12223751e-09,  -1.38472194e-08,
         6.94482484e-09,  -7.74224962e-09,   7.36062570e-07,
        -9.40351563e-06,  -2.23014657e-09,   6.74910119e-09,
        -4.93824365e-09,   8.50272392e-09,  -6.91592964e-09,
         8.88525527e-09 };

static const FLOAT coefs_mbeefvdw[25] = {
         1.17114923e+00,   1.15594371e-01,  -5.32167416e-02,
        -2.01131648e-02,   1.41417107e-03,  -6.76157938e-02,
         4.53837246e-02,  -2.22650139e-02,   1.92374554e-02,
         9.19317034e-07,   1.48659502e-02,   3.18024096e-02,
        -5.21818079e-03,   1.33707403e-07,  -5.00749348e-07,
         1.40794142e-03,  -6.08338264e-03,  -6.57949254e-07,
        -5.49909413e-08,   5.74317889e-08,   1.41530486e-04,
        -1.00478906e-07,   2.01895739e-07,   3.97324768e-09,
        -3.40722258e-09 };

typedef struct{
  int legorder;
  const FLOAT *coefs;
} mgga_x_mbeef_params;

static void
mgga_x_mbeef_init(XC(func_type) *p)
{
  mgga_x_mbeef_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_mbeef_params));
  params = (mgga_x_mbeef_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_MBEEF:
    params->legorder = 8;
    params->coefs = coefs_mbeef;
    break;
  case XC_MGGA_X_MBEEFVDW:
    params->legorder = 5;
    params->coefs = coefs_mbeefvdw;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_mbeef\n");
    exit(1);
  }
}

static inline void
transformations(int order,FLOAT p, FLOAT a,
                FLOAT *xi, FLOAT *xj, FLOAT *dxidp, FLOAT *dxjda)
{
  FLOAT tmp0, num, den, dnum, dden;
  FLOAT a2, a3;
  const FLOAT k = 6.5124; // PBEsol transformation

  /* reduced density gradient in transformation t1(s) */
  tmp0   = k + p;
  *xi    = 2.0 * p / tmp0 - 1.0;

  /* useful powers */
  a2 =  a*a; a3 = a2*a;

  /* alpha in transformation t2(a) */
  num = 1.0 - a2;
  num = -num*num*num;
  den = 1.0 + a3*(1.0 + a3);
  *xj = num / den;

  if(order < 1) return;

  *dxidp = 2.0 * k / (tmp0*tmp0);

  dnum = -a*(-6.0 + a2*(12.0 - a2*6.0));
  dden = a2*(3.0 + 6.0*a3);
  *dxjda = DFRACTION(num, dnum, den, dden);
}

static inline void
legf(int order, int legorder, const FLOAT *coefs, FLOAT xi, FLOAT xj, 
     FLOAT *f, FLOAT *dfdxi, FLOAT *dfdxj)
{
  int i, j, k, m;
  FLOAT Li[legorder], dLi[legorder], Lj[legorder], dLj[legorder];

  /* product exchange enhancement factor and derivatives */

  /* initializing */
  Li[0] = 1.0; Li[1] = xi;
  Lj[0] = 1.0; Lj[1] = xj;

  /* recursively building polynomia and their derivatives */
  for(i = 2; i < legorder; i++){
    Li[i] = 2.0 * xi * Li[i-1] - Li[i-2] - (xi * Li[i-1] - Li[i-2])/i;
    Lj[i] = 2.0 * xj * Lj[i-1] - Lj[i-2] - (xj * Lj[i-1] - Lj[i-2])/i;
  }

  /* building enhancement factor and derivatives */
  m  = 0;
  *f = 0.0;
  for(j = 0; j < legorder; j++){
    for(k = 0; k < legorder; k++){
      *f += coefs[m]*Li[k]*Lj[j];
      m++;
    }
  }

  if(order < 1) return;

  dLi[0] = 0.0; dLi[1] = 1.0;
  dLj[0] = 0.0; dLj[1] = 1.0;
  
  for(i = 2; i < legorder; i++){
    dLi[i] = i * Li[i-1] + xi * dLi[i-1];
    dLj[i] = i * Lj[i-1] + xj * dLj[i-1];    
  }
  
  m = 0;
  *dfdxi = 0.0;
  *dfdxj = 0.0;
  for(j = 0; j < legorder; j++){
    for(k = 0; k < legorder; k++){
      *dfdxi += coefs[m]*dLi[k]*Lj[j];
      *dfdxj += coefs[m]*dLj[j]*Li[k];
      m++;
    }
  }  

}

static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  FLOAT p, dpdx;
  FLOAT a, dadx, dadt;
  FLOAT xi, dxidp, xj, dxjda;
  FLOAT fx, dfdxi, dfdxj;
  mgga_x_mbeef_params *params;

  params = (mgga_x_mbeef_params *)pt->params;
  assert(params != NULL);

  p = X2S*X2S*r->x*r->x;
  a = (r->t - r->x*r->x/8.0)/K_FACTOR_C;

  transformations(r->order, p, a, &xi, &xj, &dxidp, &dxjda);
  legf(r->order, params->legorder, params->coefs, xi, xj, &fx, &dfdxi, &dfdxj);

  r->f = fx;

  if(r->order < 1) return;

  dpdx = 2.0*X2S*X2S*r->x;
  dadx = -2.0*r->x/(8.0*K_FACTOR_C);
  dadt = 1.0/K_FACTOR_C;

  r->dfdx = dfdxi*dxidp*dpdx + dfdxj*dxjda*dadx;
  r->dfdt = dfdxj*dxjda*dadt;
  r->dfdu = 0.0;
}


#include "work_mgga_x.c"


XC(func_info_type) XC(func_info_mgga_x_mbeef) = {
  XC_MGGA_X_MBEEF,
  XC_EXCHANGE,
  "mBEEF exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Wellendorff2014_144107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_mbeef_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

XC(func_info_type) XC(func_info_mgga_x_mbeefvdw) = {
  XC_MGGA_X_MBEEFVDW,
  XC_EXCHANGE,
  "mBEEF-vdW exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Lundgaard, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_mbeef_init,
  NULL, NULL, NULL,
  work_mgga_x,
};
