/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_TPSS          231 /* Tao, Perdew, Staroverov & Scuseria correlation */
#define XC_MGGA_C_TM            251 /* Tao and Mo 2016 correlation */

typedef struct{
  double beta, d;
  double C0_c[4];
} mgga_c_tpss_params;

static const mgga_c_tpss_params par_tm = {0.06672455060314922, 2.8, {0.0, 0.1, 0.32, 0.0}};

static void 
mgga_c_tpss_init(xc_func_type *p)
{
  mgga_c_tpss_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_tpss_params));
  params = (mgga_c_tpss_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_C_TPSS:
    /* default set by set_ext_params */
    break;
  case XC_MGGA_C_TM:
    memcpy(params, &par_tm, sizeof(mgga_c_tpss_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_tpss\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_beta", 0.06672455060314922, "beta"},
  {"_d", 2.8, "d"},
  {"_C0_c0", 0.53, "C0_c[0]"},
  {"_C0_c1", 0.87, "C0_c[1]"},
  {"_C0_c2", 0.50, "C0_c[2]"},
  {"_C0_c3", 2.26, "C0_c[3]"}
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_c_tpss_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_c_tpss_params *) (p->params);

  params->beta    = get_ext_param(p->info->ext_params, ext_params, 0);
  params->d       = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C0_c[0] = get_ext_param(p->info->ext_params, ext_params, 2);
  params->C0_c[1] = get_ext_param(p->info->ext_params, ext_params, 3);
  params->C0_c[2] = get_ext_param(p->info->ext_params, ext_params, 4);
  params->C0_c[3] = get_ext_param(p->info->ext_params, ext_params, 5);
}

#include "maple2c/mgga_exc/mgga_c_tpss.c"
#include "work_mgga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_tpss = {
  XC_MGGA_C_TPSS,
  XC_CORRELATION,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2003_146401, &xc_ref_Perdew2004_6898, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-23, /* densities smaller than 1e-26 give NaNs */
  6, ext_params, set_ext_params,
  mgga_c_tpss_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_tm = {
  XC_MGGA_C_TM,
  XC_CORRELATION,
  "Tao and Mo 2016 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2016_073001, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-23, /* densities smaller than 1e-26 give NaNs */
  0, NULL, NULL,
  mgga_c_tpss_init, NULL,
  NULL, NULL, work_mgga,
};
