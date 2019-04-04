/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
 Implements Tao, Perdew, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

#define XC_MGGA_X_TPSS          202 /* Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_MODTPSS       245 /* Modified Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_REVTPSS       212 /* revised Tao, Perdew, Staroverov & Scuseria exchange */
#define XC_MGGA_X_BLOC          244 /* functional with balanced localization */

typedef struct{
  double b, c, e, kappa, mu;
  double BLOC_a, BLOC_b;
} mgga_x_tpss_params;

static const mgga_x_tpss_params par_modtpss = {0.40, 1.38496, 1.37, 0.804, 0.252, 2.0, 0.0};
static const mgga_x_tpss_params par_revtpss = {0.40, 2.35203946, 2.16769874, 0.804, 0.14, 3.0, 0.0};
static const mgga_x_tpss_params par_bloc    = {0.40, 1.59096, 1.537, 0.804, 0.21951, 4.0, -3.3};

static void 
mgga_x_tpss_init(xc_func_type *p)
{
  mgga_x_tpss_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_tpss_params));
  params = (mgga_x_tpss_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_TPSS:
    /* default set by set_ext_params */
    break;
  case XC_MGGA_X_MODTPSS:
    memcpy(params, &par_modtpss, sizeof(mgga_x_tpss_params));
    break;
  case XC_MGGA_X_REVTPSS:
    memcpy(params, &par_revtpss, sizeof(mgga_x_tpss_params));
    break;
  case XC_MGGA_X_BLOC:
    memcpy(params, &par_bloc, sizeof(mgga_x_tpss_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_tpss\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_b",      0.40,    "b"},
  {"_c",      1.59096, "c"},
  {"_e",      1.537,   "e"},
  {"_kappa",  0.8040,  "Asymptotic value of the enhancement function"},
  {"_mu",     0.21951, "Coefficient of the 2nd order expansion"},
  {"_BLOC_a", 2.0,     "BLOC_a"},
  {"_BLOC_b", 0.0,     "BLOC_b"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_tpss_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_tpss_params *) (p->params);

  params->b      = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c      = get_ext_param(p->info->ext_params, ext_params, 1);
  params->e      = get_ext_param(p->info->ext_params, ext_params, 2);
  params->kappa  = get_ext_param(p->info->ext_params, ext_params, 3);
  params->mu     = get_ext_param(p->info->ext_params, ext_params, 4);
  params->BLOC_a = get_ext_param(p->info->ext_params, ext_params, 5);
  params->BLOC_b = get_ext_param(p->info->ext_params, ext_params, 6);
}

#include "maple2c/mgga_exc/mgga_x_tpss.c"
#include "work_mgga_new.c"

const xc_func_info_type xc_func_info_mgga_x_tpss = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2003_146401, &xc_ref_Perdew2004_6898, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-23,
  7, ext_params, set_ext_params,
  mgga_x_tpss_init, NULL, 
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_mgga_x_modtpss = {
  XC_MGGA_X_MODTPSS,
  XC_EXCHANGE,
  "Modified Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2007_042506, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-23,
  0, NULL, NULL,
  mgga_x_tpss_init, NULL, 
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_mgga_x_revtpss = {
  XC_MGGA_X_REVTPSS,
  XC_EXCHANGE,
  "revised Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2009_026403, &xc_ref_Perdew2009_026403_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-23,
  0, NULL, NULL,
  mgga_x_tpss_init, NULL, 
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_mgga_x_bloc = {
  XC_MGGA_X_BLOC,
  XC_EXCHANGE,
  "functional with balanced localization",
  XC_FAMILY_MGGA,
  {&xc_ref_Constantin2013_2256, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-23,
  0, NULL, NULL,
  mgga_x_tpss_init, NULL, 
  NULL, NULL, work_mgga,
};
