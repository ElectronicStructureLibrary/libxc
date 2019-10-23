/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_REVTPSS       241 /* revised TPSS correlation */
#define XC_MGGA_C_REVTM         694 /* revised Tao and Mo 2016 correlation */

typedef struct{
  double d;
  double C0_c[4];
} mgga_c_revtpss_params;

static const mgga_c_revtpss_params par_revtm = {2.8, {0.0, 0.1, 0.32, 0.0}};

static void 
mgga_c_revtpss_init(xc_func_type *p)
{
  mgga_c_revtpss_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_revtpss_params));
  params = (mgga_c_revtpss_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_C_REVTPSS:
    /* default set by set_ext_params */
    break;
  case XC_MGGA_C_REVTM:
    memcpy(params, &par_revtm, sizeof(mgga_c_revtpss_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_revtpss\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_d", 2.8, "d"},
  {"_C0_c0", 0.59,   "C0_c[0]"},
  {"_C0_c1", 0.9269, "C0_c[1]"},
  {"_C0_c2", 0.6225, "C0_c[2]"},
  {"_C0_c3", 2.1540, "C0_c[3]"}
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_c_revtpss_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_c_revtpss_params *) (p->params);

  params->d       = get_ext_param(p->info->ext_params, ext_params, 0);
  params->C0_c[0] = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C0_c[1] = get_ext_param(p->info->ext_params, ext_params, 2);
  params->C0_c[2] = get_ext_param(p->info->ext_params, ext_params, 3);
  params->C0_c[3] = get_ext_param(p->info->ext_params, ext_params, 4);
}

#include "maple2c/mgga_exc/mgga_c_revtpss.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revtpss = {
  XC_MGGA_C_REVTPSS,
  XC_CORRELATION,
  "revised TPSS correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2009_026403, &xc_ref_Perdew2009_026403_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-13, /* densities smaller than 1e-26 give NaNs */
  5, ext_params, set_ext_params,
  mgga_c_revtpss_init, NULL,
  NULL, NULL, work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_revtm = {
  XC_MGGA_C_REVTM,
  XC_CORRELATION,
  "revised Tao and Mo 2016 exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Jana2019_6356, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-13, /* densities smaller than 1e-26 give NaNs */
  0, NULL, NULL,
  mgga_c_revtpss_init, NULL,
  NULL, NULL, work_mgga
};

