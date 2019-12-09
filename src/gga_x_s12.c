/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_S12G         495 /* Swart 2012 GGA exchange                                    */
#define XC_HYB_GGA_X_S12H     496 /* Swart 2012 GGA hybrid exchange                             */

typedef struct {
  double A;
  double B;
  double C;
  double D;
  double E;
  double bx;
} gga_x_s12_params;

static void
gga_x_s12_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_s12_params));
}

static const func_params_type ext_params_s12g[] = {
  {"_A", 1.03842032, "A parameter"},
  {"_B", 1.757-1.03842032, "B parameter"},
  {"_C", 0.00403198, "C parameter"},
  {"_D", 0.00104596, "D parameter"},
  {"_E", 0.00594635, "E parameter"}
};

static void 
set_ext_params_s12g(xc_func_type *p, const double *ext_params)
{
  gga_x_s12_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_s12_params *) (p->params);

  params->A   = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B   = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C   = get_ext_param(p->info->ext_params, ext_params, 2);
  params->D   = get_ext_param(p->info->ext_params, ext_params, 3);
  params->E   = get_ext_param(p->info->ext_params, ext_params, 4);
  params->bx  = 1.0;
}

static const func_params_type ext_params_s12h[] = {
  {"_A", 1.02543951, "A parameter"},
  {"_B", 1.757-1.02543951, "B parameter"},
  {"_C", 0.00761554, "C parameter"},
  {"_D", 0.00211063, "D parameter"},
  {"_E", 0.00604672, "E parameter"},
  {"_alpha", 0.25, "Fraction of exact exchange"}
};

static void 
set_ext_params_s12h(xc_func_type *p, const double *ext_params)
{
  gga_x_s12_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_s12_params *) (p->params);

  params->A    = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B    = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C    = get_ext_param(p->info->ext_params, ext_params, 2);
  params->D    = get_ext_param(p->info->ext_params, ext_params, 3);
  params->E    = get_ext_param(p->info->ext_params, ext_params, 4);
  p->cam_alpha = get_ext_param(p->info->ext_params, ext_params, 5);
  params->bx   = 1.0 - p->cam_alpha;
  p->cam_beta  = 0.0;
  p->cam_omega = 0.0;
}

#include "maple2c/gga_exc/gga_x_s12.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_s12g = {
  XC_GGA_X_S12G,
  XC_EXCHANGE,
  "Swart 2012 GGA exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  5, ext_params_s12g, set_ext_params_s12g,
  gga_x_s12_init, NULL, 
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_s12h = {
  XC_HYB_GGA_X_S12H,
  XC_EXCHANGE,
  "Swart 2012 hybrid exchange",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  6, ext_params_s12h, set_ext_params_s12h,
  gga_x_s12_init, NULL, 
  NULL, work_gga, NULL
};
