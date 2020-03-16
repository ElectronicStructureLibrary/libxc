/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_LCGAU 708      /* Long-range Gaussian */
#define XC_GGA_X_LCGAU_CORE 709 /* Long-range Gaussian fitted to core excitations */
#define XC_GGA_X_LC2GAU 710     /* Long-range Gaussian 2 */

typedef struct{
  double a1, k1, a2, k2;
} gga_x_lcgau_params;

static void 
gga_x_lgau_init(xc_func_type *p)
{
  gga_x_lcgau_params *params;
  
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_lcgau_params));
  params = (gga_x_lcgau_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_LCGAU:
    /* default values set by set_ext_params */
    break;
  case XC_GGA_X_LCGAU_CORE:
    p->cam_omega = 0.42;
    params->a1   = 0.0335;
    params->k1   = -5.9;
    params->a2   = 0.0335;
    params->k2   = 0.0;
    break;
  case XC_GGA_X_LC2GAU:
    p->cam_omega = 0.42;
    params->a1   = 0.012;
    params->k1   = -100.0;
    params->a2   = 0.01046;
    params->k2   = 101.0;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_lcgau\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_omega", 0.42, "Screening parameter"},
  {"_a1", 0.011, "1/a multiplies the exponent of the exponential"},
  {"_k1", -18.0, "prefactor"},
  {"_a2", 0.011, "1/a multiplies the exponent of the exponential"},
  {"_k2", 0.0, "prefactor"}
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_lcgau_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_lcgau_params *) (p->params);

  p->cam_omega = get_ext_param(p->info->ext_params, ext_params, 0);
  params->a1 = get_ext_param(p->info->ext_params, ext_params, 1);
  params->k1 = get_ext_param(p->info->ext_params, ext_params, 2);
  params->a2 = get_ext_param(p->info->ext_params, ext_params, 3);
  params->k2 = get_ext_param(p->info->ext_params, ext_params, 4);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_lcgau.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lcgau = {
  XC_GGA_X_LCGAU,
  XC_EXCHANGE,
  "Long-range Gaussian",
  XC_FAMILY_GGA,
  {&xc_ref_Song2007_154109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-8,
  3, ext_params, set_ext_params,
  gga_x_lgau_init, NULL, 
  NULL, work_gga, NULL
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lcgau_core = {
  XC_GGA_X_LCGAU_CORE,
  XC_EXCHANGE,
  "Long-range Gaussian fitted to core excitations",
  XC_FAMILY_GGA,
  {&xc_ref_Song2008_184113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-8,
  0, NULL, NULL,
  gga_x_lgau_init, NULL, 
  NULL, work_gga, NULL
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lc2gau = {
  XC_GGA_X_LC2GAU,
  XC_EXCHANGE,
  "Long-range Gaussian 2",
  XC_FAMILY_GGA,
  {&xc_ref_Song2009_144108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-8,
  0, NULL, NULL,
  gga_x_lgau_init, NULL, 
  NULL, work_gga, NULL
};
