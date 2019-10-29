/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MS0          221 /* MS exchange of Sun, Xiao, and Ruzsinszky */
#define XC_MGGA_X_MS1          222 /* MS1 exchange of Sun, et al */
#define XC_MGGA_X_MS2          223 /* MS2 exchange of Sun, et al */
#define XC_HYB_MGGA_X_MS2H     224 /* MS2 hybrid exchange of Sun, et al */

typedef struct{
  double kappa, c, b;
} mgga_x_ms_params;

static void
mgga_x_ms_init(xc_func_type *p)
{
  mgga_x_ms_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_ms_params));
  params = (mgga_x_ms_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_MS0:
    /* set by set_ext_params */
    break;
  case XC_MGGA_X_MS1:
    /* set by set_ext_params */
    break;
  case XC_MGGA_X_MS2:
    /* set by set_ext_params */
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_ms\n");
    exit(1);
  }
}

static const func_params_type ms0_params[] = {
  {"_kappa",  0.29, "kappa parameter"},
  {"_c",  0.28771, "c parameter"},
  {"_b",  1.0, "exponent b"},
};

static const func_params_type ms1_params[] = {
  {"_kappa",  0.404, "kappa parameter"},
  {"_c",  0.18150, "c parameter"},
  {"_b",  1.0, "exponent b"},
};

static const func_params_type ms2_params[] = {
  {"_kappa", 0.504, "kappa parameter"},
  {"_c",  0.14601, "c parameter"},
  {"_b",  4.0, "exponent b"},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_ms_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_ms_params *) (p->params);

  params->kappa = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c     = get_ext_param(p->info->ext_params, ext_params, 1);
  params->b     = get_ext_param(p->info->ext_params, ext_params, 2);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_ms.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms0 = {
  XC_MGGA_X_MS0,
  XC_EXCHANGE,
  "MS exchange of Sun, Xiao, and Ruzsinszky",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2012_051101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-23,
  3, ms0_params, set_ext_params,
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms1 = {
  XC_MGGA_X_MS1,
  XC_EXCHANGE,
  "MS1 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-23,
  3, ms1_params, set_ext_params,
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ms2 = {
  XC_MGGA_X_MS2,
  XC_EXCHANGE,
  "MS2 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-23,
  3, ms2_params, set_ext_params,
  mgga_x_ms_init, NULL,
  NULL, NULL, work_mgga,
};

static void
hyb_mgga_x_ms2h_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_MS2};
  static double funcs_coef[1] = {0.91};

  xc_mix_init(p, 1, funcs_id, funcs_coef);
  p->cam_alpha = 0.09;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_ms2h = {
  XC_HYB_MGGA_X_MS2H,
  XC_EXCHANGE,
  "MS2 hybrid exchange of Sun, et al",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_x_ms2h_init, NULL,
  NULL, NULL, NULL
};
