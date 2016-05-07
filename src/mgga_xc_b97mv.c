/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

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

#define XC_MGGA_XC_B97M_V        254 /* Mardirossian and Head-Gordon */
#define XC_HYB_MGGA_XC_WB97M_V   531 /* Mardirossian and Head-Gordon */

#define WMAX 7
#define UMAX 5
#define NPAR 6

/* Sparse datatype for functionals */
typedef struct {
  /* Value of matrix element */
  FLOAT val;
  /* row */
  int i;
  /* column */
  int j;
} b97_parameter_t;

static const b97_parameter_t b97mv_params[][NPAR] = {
  { /* x */
    {  1.000, 0, 0},
    {  1.308, 0, 1},
    {  1.901, 0, 2},
    {  0.416, 1, 0},
    {  3.070, 1, 1},
    {  0.0,   0, 0} /* dummy entry */
  },
  { /* css */
    {  1.000, 0, 0},
    { -1.855, 0, 2},
    { -5.668, 1, 0},
    {-20.497, 3, 2},
    {-20.364, 4, 2},
    {  0.0,   0, 0} /* dummy entry */
  },
  { /* cos */
    {  1.000, 0, 0},
    {  1.573, 0, 1},
    { -6.298, 0, 3},
    {  2.535, 1, 0},
    { -6.427, 3, 2},
    {  0.0,   0, 0} /* dummy entry */
  }
};

static const b97_parameter_t wb97mv_params[][NPAR] = {
  { /* x */
    {  0.85,  0, 0},
    {  1.007, 0, 1},
    {  0.259, 1, 0},
    {  0.0,   0, 0}, /* dummy entry */
    {  0.0,   0, 0}, /* dummy entry */
    {  0.0,   0, 0}  /* dummy entry */
  },
  { /* css */
    {  0.443,  0, 0},
    { -1.437,  0, 4},
    { -4.535,  1, 0},
    { -3.39,   2, 0},
    {  4.278,  4, 3},
    {  0.0,    0, 0}  /* dummy entry */
  },
  { /* cos */
    {  1.000,  0, 0},
    {  1.358,  1, 0},
    {  2.924,  2, 0},
    { -8.812,  2, 1},
    { -1.39,   6, 0},
    {  9.142,  6, 1}
  }
};

typedef struct{
  const b97_parameter_t (*cc)[NPAR];
} mgga_xc_b97mv_params;

static void
mgga_xc_b97mv_init(XC(func_type) *p)
{
  mgga_xc_b97mv_params *params;

  assert(p != NULL);

  p->n_func_aux  = 2;
  p->func_aux    = (XC(func_type) **) malloc(2*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));
  p->func_aux[1] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_X,        XC_POLARIZED);
  XC(func_init)(p->func_aux[1], XC_LDA_C_PW_MOD, XC_POLARIZED);

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_xc_b97mv_params));
  params = (mgga_xc_b97mv_params *) p->params;

  /* Non-local correlation parameters */
  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
  
  switch(p->info->number){
  case XC_MGGA_XC_B97M_V:
    params -> cc = b97mv_params;
    break;
  case XC_HYB_MGGA_XC_WB97M_V:
    params -> cc = wb97mv_params;
    p->cam_omega =  0.3;
    XC(lda_x_set_params)(p->func_aux[0], 4.0/3.0, XC_NON_RELATIVISTIC, 0.3);
    p->cam_alpha =  1.0;
    p->cam_beta  = -(1.0 - 0.15);
    p->nlc_b = 6.0;
    p->nlc_C = 0.01;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_xc_b97mv\n");
    exit(1);
  }
}

/* tau (t) is defined here without the factor 1/2 */
static void
XC(mgga_b97mv_g)(const b97_parameter_t C[NPAR], FLOAT gamma, int order, FLOAT x, FLOAT t, FLOAT *g, FLOAT *dgdx, FLOAT *dgdt, int os)
{
  FLOAT x2, denom, u, w, dgdu, dudx, dgdw, dwdt;
  FLOAT up[UMAX+1];
  FLOAT wp[WMAX+1];
  
  int i, j;

  x2 = x*x;
  denom = 1.0 + gamma*x2;
  u = gamma*x2/denom;

  /* different definition of w for opposite spin */
  if(os)
    w = (t-1.0)/(t+1.0);
  else
    w = (2.0*K_FACTOR_C - t)/(2.0*K_FACTOR_C + t);
  
  /* calculate u^i.
     Array is offset by 1 to avoid if clause for derivative */
  up[0]=0.0;
  up[1]=1.0;
  for(i=1;i<UMAX;i++)
    up[i+1]=up[i]*u;

  /* calculate w^i.
     Array is offset by 1 to avoid if clause for derivative */
  wp[0]=0.0;
  wp[1]=1.0;
  for(i=1;i<WMAX;i++)
    wp[i+1]=wp[i]*w;
  
  /* evaluate enhancement factor */
  *g = 0.0;
  for(i=0;i<NPAR;i++)
    *g += C[i].val * wp[C[i].i+1] * up[C[i].j+1];

  if(order < 1) return;

  /* evaluate derivatives */
  dgdu = 0.0;
  for(i=0;i<NPAR;i++)
      dgdu += C[i].j * C[i].val * wp[C[i].i+1] * up[C[i].j];

  dgdw = 0.0;
  for(i=0;i<NPAR;i++)
    dgdw += C[i].i * C[i].val * wp[C[i].i] * up[C[i].j+1];
  
  dudx = 2.0*gamma*x/(denom*denom);
  /* different definition of w for opposite spin */
  if(os)
    dwdt = 2.0/((1.0 + t)*(1.0 + t));
  else
    dwdt = -4.0*K_FACTOR_C/((2.0*K_FACTOR_C + t)*(2.0*K_FACTOR_C + t));
  
  *dgdx = dgdu*dudx;
  *dgdt = dgdw*dwdt;
}

static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static const FLOAT sign[2] = {1.0, -1.0};
  static const FLOAT gamma[3] = {0.004, 0.2, 0.006}; /* Xs, Css, Cab */
  static const FLOAT tmin = 1e-10;
  /* run opposite spin correlation? */
  int cos=1;
  
  mgga_xc_b97mv_params *params;
  
  XC(lda_work_t) lda_x[3], lda_pw[3];
  FLOAT opz;
  FLOAT gx, dgxdx, dgxdt;
  FLOAT gcss, dgcssdx, dgcssdt;
  FLOAT gcos, dgcosdx, dgcosdt;
  FLOAT aux, aux12, x_avg, dx_avgdxs[2];
  FLOAT t[2], t_avg;
  int is;

  /* get the parameters */
  params = (mgga_xc_b97mv_params *) pt->params;
  
  /* first we get the parallel and perpendicular LDAs */
  XC(lda_stoll) (pt->func_aux[0], XC(lda_x_func),    r->dens, r->zeta, r->order, lda_x);
  XC(lda_stoll) (pt->func_aux[1], XC(lda_c_pw_func), r->dens, r->zeta, r->order, lda_pw);

  /* initialize to zero */
  r->f = 0.0;
  if(r->order >= 1){
    r->dfdrs = r->dfdz = r->dfdxs[0] = r->dfdxs[1] = r->dfdxt = 0.0;
    r->dfdus[0] = r->dfdus[1] = r->dfdts[0] = r->dfdts[1] = 0.0;
  }
  if(r->order >= 2){
    r->d2fdrs2 = r->d2fdrsz = r->d2fdrsxt = r->d2fdrsxs[0] = r->d2fdrsxs[1] = 0.0;
    r->d2fdz2 = r->d2fdzxt = r->d2fdzxs[0] = r->d2fdzxs[1] = r->d2fdxt2 = 0.0;
    r->d2fdxtxs[0] = r->d2fdxtxs[1] = r->d2fdxs2[0] = r->d2fdxs2[1] = r->d2fdxs2[2] = 0.0;
  }

  /* now we calculate the g functions for exchange and parallel correlation */
  for(is = 0; is < 2; is++){
    opz   = 1.0 + sign[is]*r->zeta;

    if(r->dens*opz < 2.0*pt->info->min_dens) {
      cos=0;
      continue;
    }
    
    XC(mgga_b97mv_g)(params->cc[0], gamma[0], r->order, r->xs[is], 2.0*r->ts[is], &gx, &dgxdx, &dgxdt, 0);
    XC(mgga_b97mv_g)(params->cc[1], gamma[1], r->order, r->xs[is], 2.0*r->ts[is], &gcss, &dgcssdx, &dgcssdt, 0);
    r->f += lda_x[is].zk*gx + lda_pw[is].zk*gcss;

    if(r->order < 1) continue;

    r->dfdrs     += lda_x[is].dedrs *  gx   + lda_pw[is].dedrs *  gcss;
    r->dfdz      += lda_x[is].dedz  *  gx   + lda_pw[is].dedz  *  gcss;
    r->dfdxs[is] += lda_x[is].zk    * dgxdx + lda_pw[is].zk    * dgcssdx;
    r->dfdts[is] += lda_x[is].zk * (2.0*dgxdt) + lda_pw[is].zk * (2.0*dgcssdt);
  }

  /* and now we add the opposite-spin contribution */
  if(cos) {
    aux = r->xs[0]*r->xs[0] + r->xs[1]*r->xs[1];
    aux12 = SQRT(aux);
    x_avg = aux12/M_SQRT2;

    for(is=0;is<2;is++)
      t[is] =  (r->ts[is] > tmin) ? K_FACTOR_C/r->ts[is] : 0.0;
    t_avg = (t[0] + t[1])/2.0;
    
    XC(mgga_b97mv_g)(params->cc[2], gamma[2], r->order, x_avg, t_avg, &gcos, &dgcosdx, &dgcosdt, 1);
    r->f += lda_pw[2].zk*gcos;
    
    if(r->order < 1) return;
    
    dx_avgdxs[0] = r->xs[0]/(aux12*M_SQRT2);
    dx_avgdxs[1] = r->xs[1]/(aux12*M_SQRT2);
    
    r->dfdrs    += lda_pw[2].dedrs *  gcos;
    r->dfdz     += lda_pw[2].dedz  *  gcos;
    r->dfdxs[0] += lda_pw[2].zk    * dgcosdx * dx_avgdxs[0];
    r->dfdxs[1] += lda_pw[2].zk    * dgcosdx * dx_avgdxs[1];
    r->dfdts[0] += lda_pw[2].zk    * dgcosdt * (-t[0]*t[0]/(2.0*K_FACTOR_C));
    r->dfdts[1] += lda_pw[2].zk    * dgcosdt * (-t[1]*t[1]/(2.0*K_FACTOR_C));
  }
}

#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_xc_b97m_v) = {
  XC_MGGA_XC_B97M_V,
  XC_EXCHANGE_CORRELATION,
  "B97M-V exchange-correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Mardirossian2015_074111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_xc_b97mv_init,
  NULL, NULL, NULL,
  work_mgga_c,
};

const XC(func_info_type) XC(func_info_hyb_mgga_xc_wb97m_v) = {
  XC_HYB_MGGA_XC_WB97M_V,
  XC_EXCHANGE_CORRELATION,
  "wB97M-V exchange-correlation functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Mardirossian2016, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_xc_b97mv_init,
  NULL, NULL, NULL,
  work_mgga_c,
};
