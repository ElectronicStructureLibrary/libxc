/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_KT1          145 /* Exchange part of Keal and Tozer version 1 */
#define XC_GGA_XC_KT1         167 /* Keal and Tozer version 1                  */
#define XC_GGA_XC_KT2         146 /* Keal and Tozer version 2                  */
#define XC_GGA_XC_KT3         587 /* Keal and Tozer version 3                  */

typedef struct{
  double gamma, delta;
} gga_x_kt_params;

static void 
gga_x_kt_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_kt_params));
}

static const func_params_type ext_params[] = {
  {"_gamma", -0.006, "gamma"},
  {"_delta",    0.1, "delta"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_kt_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_kt_params *) (p->params);

  params->gamma = get_ext_param(p->info->ext_params, ext_params, 0);
  params->delta = get_ext_param(p->info->ext_params, ext_params, 1);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_kt.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_kt1 = {
  XC_GGA_X_KT1,
  XC_EXCHANGE,
  "Exchange part of Keal and Tozer version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  2, ext_params, set_ext_params,
  gga_x_kt_init, NULL, 
  NULL, work_gga, NULL
};


static void
gga_xc_kt1_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);  
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt1 = {
  XC_GGA_XC_KT1,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  0, NULL, NULL,
  gga_xc_kt1_init, NULL, 
  NULL, NULL, NULL
};


static void
gga_xc_kt2_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[3] = {1.07173 - 1.0, 1.0, 0.576727};

  xc_mix_init(p, 3, funcs_id, funcs_coef);  
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt2 = {
  XC_GGA_XC_KT2,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 2",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  0, NULL, NULL,
  gga_xc_kt2_init, NULL, 
  NULL, NULL, NULL
};


static void
gga_xc_kt3_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_C_LYP, XC_GGA_X_KT1, XC_GGA_X_OPTX};
  double funcs_coef[4];

  /* Equation (2) */
  static const double
    alpha = 1.092,
    beta  = 0.864409,
    par_kt[2] = {-0.004, 0.1},
    eps   = -0.925452,
    a1    = 1.05151, /* these are the OPTX coefficients */
    a2    = 1.43169;

  funcs_coef[0] = alpha - eps*a1/a2;
  funcs_coef[1] = beta;
  funcs_coef[2] = 1.0;
  funcs_coef[3] = eps/a2;

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  set_ext_params(p->func_aux[2], par_kt);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_kt3 = {
  XC_GGA_XC_KT3,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 3",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2004_5654, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  0, NULL, NULL,
  gga_xc_kt3_init, NULL, 
  NULL, NULL, NULL
};
