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

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_PBE_SOL      116 /* Perdew, Burke & Ernzerhof exchange (solids)    */
#define XC_GGA_X_XPBE         123 /* xPBE reparametrization by Xu & Goddard         */
#define XC_GGA_X_PBE_JSJR     126 /* JSJR reparametrization by Pedroza, Silva & Capelle */
#define XC_GGA_X_PBEK1_VDW    140 /* PBE reparametrization for vdW                  */
#define XC_GGA_X_RGE2         142 /* Regularized PBE                                */
#define XC_GGA_X_APBE         184 /* mu fixed from the semiclassical neutral atom   */
#define XC_GGA_X_PBEINT        60 /* PBE for hybrid interfaces                      */
#define XC_GGA_X_PBE_TCA       59 /* PBE revised by Tognetti et al                  */
#define XC_GGA_X_PBE_MOL       49 /* Del Campo, Gazquez, Trickey and Vela (PBE-like) */
#define XC_GGA_X_LAMBDA_LO_N   45 /* lambda_LO(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_CH_N   44 /* lambda_CH(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_OC2_N  40 /* lambda_OC2(N) version of PBE                   */
#define XC_GGA_X_BCGP          38 /* Burke, Cancio, Gould, and Pittalis             */
#define XC_GGA_X_PBEFE        265 /* PBE for formation energies                     */

#define XC_GGA_K_APBE         185 /* mu fixed from the semiclassical neutral atom   */
#define XC_GGA_K_TW1          187 /* Tran and Wesolowski set 1 (Table II)           */
#define XC_GGA_K_TW2          188 /* Tran and Wesolowski set 2 (Table II)           */
#define XC_GGA_K_TW3          189 /* Tran and Wesolowski set 3 (Table II)           */
#define XC_GGA_K_TW4          190 /* Tran and Wesolowski set 4 (Table II)           */
#define XC_GGA_K_REVAPBE       55 /* revised APBE                                   */
#define XC_GGA_K_APBEINT       54 /* interpolated version of APBE                   */
#define XC_GGA_K_REVAPBEINT    53 /* interpolated version of REVAPBE                */


typedef struct{
  FLOAT kappa, mu;

  /* parameters only used for PBEint and similar functionals */
  FLOAT alpha, muPBE, muGE;

  /* parameter used in the Odashima & Capelle versions */
  FLOAT lambda;
} gga_x_pbe_params;


static void 
gga_x_pbe_init(XC(func_type) *p)
{
  static const FLOAT kappa[] = {
    0.8040,                 /* original PBE */
    1.245,                  /* PBE R       */
    0.8040,                 /* PBE sol     */
    0.91954,                /* xPBE        */
    0.8040,                 /* PBE_JSJR    */
    1.0,                    /* PBEK1_VDW   */
    0.8040,                 /* RGE2        */
    0.8040,                 /* APBE (X)    */
    0.8040,                 /* APBE (K)    */
    0.8209,                 /* TW1         */
    0.6774,                 /* TW2         */
    0.8438,                 /* TW3         */
    0.8589,                 /* TW4         */
    0.8040,                 /* PBEint      */
    1.227,                  /* PBE TCA     */
    1.245,                  /* revAPBE (K) */
    0.8040,                 /* APBEINT (K) */
    1.245,                  /* revAPBEINT (K) */
    0.8040,                 /* PBEmol    */
    2.273/M_CBRT2 - 1.0,    /* LAMBDA_LO(N)  */
    2.215/M_CBRT2 - 1.0,    /* LAMBDA_CH(N)  */
    2.00 /M_CBRT2 - 1.0,    /* LAMBDA_OC2(N) */
    0.8040,                 /* BCGP (X)    */
    0.437,                  /* PBEfe       */
  };

  static const FLOAT mu[] = {
    0.2195149727645171,     /* PBE: mu = beta*pi^2/3, beta = 0.06672455060314922 */
    0.2195149727645171,     /* PBE rev: as PBE */
    10.0/81.0,              /* PBE sol */
    0.23214,                /* xPBE */
    0.046*M_PI*M_PI/3.0,    /* PBE_JSJR */
    0.2195149727645171,     /* PBEK1_VDW: as PBE */
    10.0/81.0,              /* RGE2      */
    0.260,                  /* APBE (X)  */
    0.23889,                /* APBE (K)  */
    0.2335,                 /* TW1       */
    0.2371,                 /* TW2       */
    0.2319,                 /* TW3       */
    0.2309,                 /* TW4       */
    0.0,                    /* PBEint (to be set later) */
    0.2195149727645171,     /* PBE TCA: as PBE */
    0.23889,                /* revAPBE (K)  */
    0.0,                    /* APBEINT (K) (to be set later) */
    0.0,                    /* REVAPBEINT (K) (to be set later) */
    0.27583,                /* PBEmol        */
    0.2195149727645171,     /* LAMBDA_LO(N)  */
    0.2195149727645171,     /* LAMBDA_CH(N)  */
    0.2195149727645171,     /* LAMBDA_OC2(N) */
    0.249,                  /* BCGP (X)    */
    0.346,                  /* PBEfe       */
  };

  gga_x_pbe_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_pbe_params));
  params = (gga_x_pbe_params *) (p->params);
 
  /* initialize structure */
  params->alpha  = 0.0;
  params->muPBE  = 0.0;
  params->muGE   = 0.0;
  params->lambda = 0.0;

  switch(p->info->number){
  case XC_GGA_X_PBE:        p->func = 0;  break;
  case XC_GGA_X_PBE_R:      p->func = 1;  break;
  case XC_GGA_X_PBE_SOL:    p->func = 2;  break;
  case XC_GGA_X_XPBE:       p->func = 3;  break;
  case XC_GGA_X_PBE_JSJR:   p->func = 4;  break;
  case XC_GGA_X_PBEK1_VDW:  p->func = 5;  break;
  case XC_GGA_X_RGE2:       p->func = 6;  break;
  case XC_GGA_X_APBE:       p->func = 7;  break;
  case XC_GGA_K_APBE:       p->func = 8;  break;
  case XC_GGA_K_TW1:        p->func = 9;  break;
  case XC_GGA_K_TW2:        p->func = 10; break;
  case XC_GGA_K_TW3:        p->func = 11; break;
  case XC_GGA_K_TW4:        p->func = 12; break;
  case XC_GGA_X_PBEINT: {
    p->func  = 13;
    params->alpha = 0.197;
    params->muPBE = 0.2195149727645171;
    params->muGE  = 10.0/81.0;
    break;
  }
  case XC_GGA_X_PBE_TCA:    p->func = 14; break;
  case XC_GGA_K_REVAPBE:    p->func = 15; break;
  case XC_GGA_K_APBEINT: {
    p->func  = 16;
    params->alpha = 5.0/3.0;
    params->muPBE = 0.23899;
    params->muGE  = 5.0/27.0;
    break;
  }
  case XC_GGA_K_REVAPBEINT: { /* equal to the previous one */
    p->func  = 17;
    params->alpha = 5.0/3.0;
    params->muPBE = 0.23899;
    params->muGE  = 5.0/27.0;
    break;
  }
  case XC_GGA_X_PBE_MOL:    p->func = 18; break;

  case XC_GGA_X_LAMBDA_LO_N: {
    p->func = 19;
    params->lambda = 2.273;
    break;
  }
  case XC_GGA_X_LAMBDA_CH_N: {
    p->func = 20;
    params->lambda = 2.215;
    break;
  }
  case XC_GGA_X_LAMBDA_OC2_N: {
    p->func = 21;
    params->lambda = 2.00;
    break;
  }
  case XC_GGA_X_BCGP:    p->func = 22; break;
  case XC_GGA_X_PBEFE:   p->func = 23; break;

  default:{
    fprintf(stderr, "Internal error in gga_x_pbe\n");
    exit(1);
  }}

  XC(gga_x_pbe_set_params)(p, kappa[p->func], mu[p->func]);
}


void 
XC(gga_x_pbe_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT mu)
{
  gga_x_pbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  params->kappa = kappa;
  params->mu    = mu;
}


void 
XC(gga_x_lambda_set_params)(XC(func_type) *p, FLOAT N)
{
  const FLOAT lambda_1 = 1.48;

  gga_x_pbe_params *params;
  FLOAT lambda;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  lambda = (1.0 - 1.0/N)*params->lambda + lambda_1/N;
  params->kappa = lambda/M_CBRT2 - 1.0;
}


void XC(gga_x_pbe_enhance) 
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  gga_x_pbe_params *params;
  FLOAT kappa, auxmu, mu, dmu, d2mu, d3mu, ss, ss2, f0, df0, d2f0, d3f0;

  assert(p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  kappa = params->kappa;

  ss  = X2S*x;
  ss2 = ss*ss;
 
  if(params->alpha != 0.0){ /* PBEint and related functionals */
    auxmu = 1.0 + params->alpha*ss2;
    mu = params->muGE + (params->muPBE - params->muGE) * params->alpha*ss2/auxmu;
  }else
    mu = params->mu;

  f0 = kappa + mu*ss2;
  if(p->info->number == XC_GGA_X_RGE2)
    f0 += mu*mu*ss2*ss2/kappa;

  *f = 1.0 + kappa*(1.0 - kappa/f0);

  if(order < 1) return;

  if(params->alpha != 0.0){ /* PBEint and related functionals */
    dmu = (params->muPBE - params->muGE) * 2.0*params->alpha*ss/(auxmu*auxmu);
  }else
    dmu = 0.0;

  df0 = 2.0*mu*ss + dmu*ss2;
  if(p->info->number == XC_GGA_X_RGE2)
    df0 += 4.0*mu*mu*ss2*ss/kappa;

  *dfdx  = X2S*kappa*kappa*df0/(f0*f0);

  if(order < 2) return;

  if(params->alpha != 0.0) /* PBEint and related functionals */
    d2mu = (params->muPBE - params->muGE) * 
      2.0*params->alpha*(1.0 - 3.0*params->alpha*ss2)/(auxmu*auxmu*auxmu);
  else
    d2mu = 0.0;

  d2f0 = 2.0*mu + 4.0*dmu*ss + d2mu*ss2;
  if(p->info->number == XC_GGA_X_RGE2)
    d2f0 += 4.0*3.0*mu*mu*ss2/kappa;

  *d2fdx2 = -X2S*X2S*kappa*kappa*(2.0*df0*df0 - d2f0*f0)/(f0*f0*f0);

  if(order < 3) return;

  if(params->alpha != 0.0) /* PBEint and related functionals */
    d3mu = (params->muPBE - params->muGE) * 
      24.0*params->alpha*params->alpha*ss*(-1.0 + params->alpha*ss2)/(auxmu*auxmu*auxmu*auxmu);
  else
    d3mu = 0.0;  

  d3f0 = 6.0*dmu + 6.0*ss*d2mu + ss2*d3mu;
  if(p->info->number == XC_GGA_X_RGE2)
    d3f0 += 4.0*3.0*2.0*mu*mu*ss/kappa;

  *d3fdx3 = X2S*X2S*X2S*kappa*kappa*(6.0*df0*df0*df0 - 6.0*f0*df0*d2f0 + f0*f0*d3f0)/(f0*f0*f0*f0);
  
}


#define func XC(gga_x_pbe_enhance)
#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_pbe) = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbe_r) = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Revised PBE from Zhang & Yang",
  XC_FAMILY_GGA,
  {&xc_ref_Zhang1998_890, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbe_sol) = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew2008_136406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_xpbe) = {
  XC_GGA_X_XPBE,
  XC_EXCHANGE,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  {&xc_ref_Xu2004_4068, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbe_jsjr) = {
  XC_GGA_X_PBE_JSJR,
  XC_EXCHANGE,
  "Reparametrized PBE by Pedroza, Silva & Capelle",
  XC_FAMILY_GGA,
  {&xc_ref_Pedroza2009_201106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbek1_vdw) = {
  XC_GGA_X_PBEK1_VDW,
  XC_EXCHANGE,
  "Reparametrized PBE for vdW",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2010_022201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_rge2) = {
  XC_GGA_X_RGE2,
  XC_EXCHANGE,
  "Regularized PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Ruzsinszky2009_763, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_apbe) = {
  XC_GGA_X_APBE,
  XC_EXCHANGE,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbeint) = {
  XC_GGA_X_PBEINT,
  XC_EXCHANGE,
  "PBE for hybrid interfaces",
  XC_FAMILY_GGA,
  {&xc_ref_Fabiano2010_113104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};


const XC(func_info_type) XC(func_info_gga_x_pbe_tca) = {
  XC_GGA_X_PBE_TCA,
  XC_EXCHANGE,
  "PBE revised by Tognetti et al",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2008_536, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_lambda_lo_n) = {
  XC_GGA_X_LAMBDA_LO_N,
  XC_EXCHANGE,
  "lambda_LO(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_lambda_ch_n) = {
  XC_GGA_X_LAMBDA_CH_N,
  XC_EXCHANGE,
  "lambda_CH(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_lambda_oc2_n) = {
  XC_GGA_X_LAMBDA_OC2_N,
  XC_EXCHANGE,
  "lambda_OC2(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbe_mol) = {
  XC_GGA_X_PBE_MOL,
  XC_EXCHANGE,
  "Reparametrized PBE by del Campo, Gazquez, Trickey & Vela",
  XC_FAMILY_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_bcgp) = {
  XC_GGA_X_BCGP,
  XC_EXCHANGE,
  "Burke, Cancio, Gould, and Pittalis",
  XC_FAMILY_GGA,
  {&xc_ref_Burke2014_4834, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_pbefe) = {
  XC_GGA_X_PBEFE,
  XC_EXCHANGE,
  "PBE for formation energies",
  XC_FAMILY_GGA,
  {&xc_ref_Perez2015_3844, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_apbe) = {
  XC_GGA_K_APBE,
  XC_KINETIC,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_revapbe) = {
  XC_GGA_K_REVAPBE,
  XC_KINETIC,
  "revised APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_tw1) = {
  XC_GGA_K_TW1,
  XC_KINETIC,
  "Tran and Wesolowski set 1 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_tw2) = {
  XC_GGA_K_TW2,
  XC_KINETIC,
  "Tran and Wesolowski set 2 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_tw3) = {
  XC_GGA_K_TW3,
  XC_KINETIC,
  "Tran and Wesolowski set 3 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_tw4) = {
  XC_GGA_K_TW4,
  XC_KINETIC,
  "Tran and Wesolowski set 4 (Table II)",
  XC_FAMILY_GGA,
  {&xc_ref_Tran2002_441, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_apbeint) = {
  XC_GGA_K_APBEINT,
  XC_KINETIC,
  "interpolated version of APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_revapbeint) = {
  XC_GGA_K_REVAPBEINT,
  XC_KINETIC,
  "interpolated version of revAPBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init,
  NULL, NULL,
  work_gga_k,
  NULL
};
