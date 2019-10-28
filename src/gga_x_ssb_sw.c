/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_SSB_SW       90  /* Swart, Sola and Bickelhaupt correction to PBE  */
#define XC_GGA_X_SSB          91  /* Swart, Sola and Bickelhaupt  */
#define XC_GGA_X_SSB_D        92  /* Swart, Sola and Bickelhaupt dispersion  */

typedef struct{
  double A, B, C, D, E;
} gga_x_ssb_sw_params;


static void 
gga_x_ssb_sw_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ssb_sw_params));
}

static const func_params_type ext_params[] = {
  {"_A", 1.0515,   "Constant s limit"},
  {"_B", 0.191458, "B s^2/(1 + C s^2)"},
  {"_C", 0.254443, "B s^2/(1 + C s^2)"},
  {"_D", 0.180708, "D s^2/(1 + E s^4)"},
  {"_E", 4.036674, "D s^2/(1 + E s^4)"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_ssb_sw_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ssb_sw_params *) (p->params);

  params->A = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C = get_ext_param(p->info->ext_params, ext_params, 2);
  params->D = get_ext_param(p->info->ext_params, ext_params, 3);
  params->E = get_ext_param(p->info->ext_params, ext_params, 4);
}

#include "maple2c/gga_exc/gga_x_ssb_sw.c"
#include "work_gga.c"

const xc_func_info_type xc_func_info_gga_x_ssb_sw = {
  XC_GGA_X_SSB_SW,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt correction to PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_69, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-22,
  5, ext_params, set_ext_params,
  gga_x_ssb_sw_init, NULL, 
  NULL, work_gga, NULL
};

static void
gga_x_ssb_init(xc_func_type *p)
{
  static const double u = -1.205643, F = 0.995010, B = 0.137574;

  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_SSB_SW, XC_GGA_X_KT1};
  static double funcs_coef[3] = {-1.0, 1.0, 1.0};

  static double par_x_ssb_sw[] = {1.071769, 0.137574, 0.187883, 0.137574*(1.0 + 1.205643), 6.635315};
  static double par_x_kt[] = {-1, 0.1};
  par_x_kt[0] = u*F*X_FACTOR_C*B*(X2S*X2S);
  
  xc_mix_init(p, 3, funcs_id, funcs_coef);  

  xc_func_set_ext_params(p->func_aux[1], par_x_ssb_sw);
  xc_func_set_ext_params(p->func_aux[2], par_x_kt);
}


const xc_func_info_type xc_func_info_gga_x_ssb = {
  XC_GGA_X_SSB,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  0, NULL, NULL,
  gga_x_ssb_init, NULL, 
  NULL, NULL, NULL
};


static void
gga_x_ssb_d_init(xc_func_type *p)
{
  static const double u = -0.749940, F = 0.949488, B = 0.197465;

  static int    funcs_id  [3] = {XC_LDA_X, XC_GGA_X_SSB_SW, XC_GGA_X_KT1};
  static double funcs_coef[3] = {-1.0, 1.0, 1.0};

  static double par_x_ssb_sw[] = {1.079966, 0.197465, 0.272729, 0.197465*(1.0 + 0.749940), 5.873645};
  static double par_x_kt[] = {-1, 0.1};
  par_x_kt[0] = u*F*X_FACTOR_C*B*(X2S*X2S);
  
  xc_mix_init(p, 3, funcs_id, funcs_coef);  

  xc_func_set_ext_params(p->func_aux[1], par_x_ssb_sw);
  xc_func_set_ext_params(p->func_aux[2], par_x_kt);
}

const xc_func_info_type xc_func_info_gga_x_ssb_d = {
  XC_GGA_X_SSB_D,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt dispersion",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-23,
  0, NULL, NULL,
  gga_x_ssb_d_init, NULL, 
  NULL, NULL, NULL
};


