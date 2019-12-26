/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_PBE_SOL      116 /* Perdew, Burke & Ernzerhof exchange (solids)    */
#define XC_GGA_X_XPBE         123 /* xPBE reparametrization by Xu & Goddard         */
#define XC_GGA_X_PBE_JSJR     126 /* JSJR reparametrization by Pedroza, Silva & Capelle */
#define XC_GGA_X_PBEK1_VDW    140 /* PBE reparametrization for vdW                  */
#define XC_GGA_X_APBE         184 /* mu fixed from the semiclassical neutral atom   */
#define XC_GGA_X_PBE_TCA       59 /* PBE revised by Tognetti et al                  */
#define XC_GGA_X_PBE_MOL       49 /* Del Campo, Gazquez, Trickey and Vela (PBE-like) */
#define XC_GGA_X_LAMBDA_LO_N   45 /* lambda_LO(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_CH_N   44 /* lambda_CH(N) version of PBE                    */
#define XC_GGA_X_LAMBDA_OC2_N  40 /* lambda_OC2(N) version of PBE                   */
#define XC_GGA_X_BCGP          38 /* Burke, Cancio, Gould, and Pittalis             */
#define XC_GGA_X_PBEFE        265 /* PBE for formation energies                     */

typedef struct{
  double kappa, mu;
  double lambda;   /* parameter used in the Odashima & Capelle versions */
} gga_x_pbe_params;


static void 
gga_x_pbe_init(xc_func_type *p)
{
  gga_x_pbe_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_params));
  params = (gga_x_pbe_params *) (p->params);
 
  params->lambda = 0.0;

  switch(p->info->number){
  case XC_GGA_X_PBE:
    /* default set by set_ext_params */
    break;
  case XC_GGA_X_PBE_R:
    params->kappa = 1.245;
    params->mu    = MU_PBE;
    break;
  case XC_GGA_X_PBE_SOL:
    params->kappa = 0.804;
    params->mu    = MU_GE;
    break;
  case XC_GGA_X_XPBE:
    params->kappa = 0.91954;
    params->mu    = 0.23214;
    break;
  case XC_GGA_X_PBE_JSJR:
    params->kappa = 0.8040;
    params->mu    = 0.046*M_PI*M_PI/3.0;
    break;
  case XC_GGA_X_PBEK1_VDW:
    params->kappa = 1.0;
    params->mu    = MU_PBE;
    break;
  case XC_GGA_X_APBE:
    params->kappa = 0.8040;
    params->mu    = 0.260;
    break;
  case XC_GGA_X_PBE_TCA:
    params->kappa = 1.227;
    params->mu    = MU_PBE;
    break;
  case XC_GGA_X_PBE_MOL:
    params->kappa = 0.8040;
    params->mu    = 0.27583;
    break;
  case XC_GGA_X_LAMBDA_LO_N:
    params->kappa = -1.0;
    params->mu    = MU_PBE;
    params->lambda = 2.273;
    break;
  case XC_GGA_X_LAMBDA_CH_N:
    params->kappa = -1.0;
    params->mu    = MU_PBE;
    params->lambda = 2.215;
    break;
  case XC_GGA_X_LAMBDA_OC2_N:
    params->kappa = -1.0;
    params->mu    = MU_PBE;
    params->lambda = 2.00;
    break;
  case XC_GGA_X_BCGP:
    params->kappa = 0.8040;
    params->mu    = 0.249;
    break;
  case XC_GGA_X_PBEFE:
    params->kappa = 0.437;
    params->mu    = 0.346;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_pbe\n");
    exit(1);
  }
}

/* PBE: mu = beta*pi^2/3, beta = 0.06672455060314922 */
static const func_params_type ext_params_PBE[] = {
  {"_kappa", 0.8040, "Asymptotic value of the enhancement function"},
  {"_mu",    MU_PBE, "Coefficient of the 2nd order expansion"},
};

static void 
set_ext_params_PBE(xc_func_type *p, const double *ext_params)
{
  gga_x_pbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  params->kappa = get_ext_param(p->info->ext_params, ext_params, 0);
  params->mu    = get_ext_param(p->info->ext_params, ext_params, 1);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pbe.c"
#include "work_gga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  2, ext_params_PBE, set_ext_params_PBE,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_r = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Revised PBE from Zhang & Yang",
  XC_FAMILY_GGA,
  {&xc_ref_Zhang1998_890, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_sol = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew2008_136406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_xpbe = {
  XC_GGA_X_XPBE,
  XC_EXCHANGE,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  {&xc_ref_Xu2004_4068, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_jsjr = {
  XC_GGA_X_PBE_JSJR,
  XC_EXCHANGE,
  "Reparametrized PBE by Pedroza, Silva & Capelle",
  XC_FAMILY_GGA,
  {&xc_ref_Pedroza2009_201106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbek1_vdw = {
  XC_GGA_X_PBEK1_VDW,
  XC_EXCHANGE,
  "Reparametrized PBE for vdW",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2010_022201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_apbe = {
  XC_GGA_X_APBE,
  XC_EXCHANGE,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_tca = {
  XC_GGA_X_PBE_TCA,
  XC_EXCHANGE,
  "PBE revised by Tognetti et al",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2008_536, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

static const func_params_type ext_params_N[] = {
  {"N", 1e23, "Number of electrons"},
};

static void 
set_ext_params_N(xc_func_type *p, const double *ext_params)
{
  const double lambda_1 = 1.48;

  gga_x_pbe_params *params;
  double lambda, ff;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbe_params *) (p->params);

  ff = (ext_params == NULL) ? p->info->ext_params[0].value : ext_params[0];

  lambda = (1.0 - 1.0/ff)*params->lambda + lambda_1/ff;
  params->kappa = lambda/M_CBRT2 - 1.0;
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_lo_n = {
  XC_GGA_X_LAMBDA_LO_N,
  XC_EXCHANGE,
  "lambda_LO(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  1, ext_params_N, set_ext_params_N,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_ch_n = {
  XC_GGA_X_LAMBDA_CH_N,
  XC_EXCHANGE,
  "lambda_CH(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  1, ext_params_N, set_ext_params_N,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lambda_oc2_n = {
  XC_GGA_X_LAMBDA_OC2_N,
  XC_EXCHANGE,
  "lambda_OC2(N) version of PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Odashima2009_798, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  1, ext_params_N, set_ext_params_N,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbe_mol = {
  XC_GGA_X_PBE_MOL,
  XC_EXCHANGE,
  "Reparametrized PBE by del Campo, Gazquez, Trickey & Vela",
  XC_FAMILY_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_bcgp = {
  XC_GGA_X_BCGP,
  XC_EXCHANGE,
  "Burke, Cancio, Gould, and Pittalis",
  XC_FAMILY_GGA,
  {&xc_ref_Burke2014_4834, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pbefe = {
  XC_GGA_X_PBEFE,
  XC_EXCHANGE,
  "PBE for formation energies",
  XC_FAMILY_GGA,
  {&xc_ref_Perez2015_3844, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga, NULL
};
