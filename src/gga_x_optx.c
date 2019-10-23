/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_OPTX         110 /* Handy & Cohen OPTX 01                          */

typedef struct{
  double a, b, gamma;
} gga_x_optx_params;


static void 
gga_x_optx_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_optx_params));
}

static const func_params_type ext_params[] = {
  {"_a", 1.05151, "a"},
  {"_b", 1.43169/X_FACTOR_C, "b"},
  {"_gamma", 0.006, "gamma"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_optx_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_optx_params *) (p->params);

  params->a     = get_ext_param(p->info->ext_params, ext_params, 0);
  params->b     = get_ext_param(p->info->ext_params, ext_params, 1);
  params->gamma = get_ext_param(p->info->ext_params, ext_params, 2);
}

#include "maple2c/gga_exc/gga_x_optx.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  {&xc_ref_Handy2001_403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-22,
  3, ext_params, set_ext_params,
  gga_x_optx_init, NULL, 
  NULL, work_gga, NULL
};
