/*
 Copyright (C) 2006-2021 M.A.L. Marques
                    2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_APBE_BFIT        830 /* APBE with a fit of beta */
#define XC_GGA_C_PBEINT_BFIT      831 /* PBEINT with a fit of beta */
#define XC_HYB_GGA_XC_LC_WAPBE    832 /* LC hybrid based on APBE */
#define XC_HYB_GGA_XC_LC_WPBEINT  833 /* LC hybrid based on PBEINT */
#define XC_HYB_GGA_XC_LC_WSG4     834 /* LC hybrid based on SG4 */

typedef struct{
  double bfit_a1, bfit_a2, bfit_a3, bfit_a4, bfit_a5;
  double bfit_mux, bfit_omega;
} gga_c_apbe_bfit_params;


static void gga_c_apbe_bfit_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_apbe_bfit_params));
}

#define APBE_BFIT_N_PAR 7
static const char  *apbe_bfit_names[APBE_BFIT_N_PAR]  =
  {
   "a1", "a2", "a3", "a4", "a5",
   "_mux", "_omega"
  };
static const char  *apbe_bfit_desc[APBE_BFIT_N_PAR]   =
  {
   "Pade fit: a1", "Pade fit: a2", "Pade fit: a3", "Pade fit: a4", "Pade fit: a5",
   "mu of the corresponding exchange functional",
   "range separation parameter"
  };
static const double apbe_bfit_values[APBE_BFIT_N_PAR] =
  {0.06929609, 0.02090877, 73.63025684, 3.84513730, 0.00000049, 0.26, 0.37};
static const double pbeint_bfit_values[APBE_BFIT_N_PAR] =
  {0.00000000, 0.06413244, 27.06803466, 3.61233368, 0.00005694, 10./81., 0.52};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_apbe_bfit.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_apbe_bfit = {
  XC_GGA_C_APBE_BFIT,
  XC_CORRELATION,
  "APBE with a fit of beta",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {APBE_BFIT_N_PAR, apbe_bfit_names, apbe_bfit_desc, apbe_bfit_values, set_ext_params_cpy},
  gga_c_apbe_bfit_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_pbeint_bfit = {
  XC_GGA_C_PBEINT_BFIT,
  XC_CORRELATION,
  "PBEINT with a fit of beta",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {APBE_BFIT_N_PAR, apbe_bfit_names, apbe_bfit_desc, pbeint_bfit_values, set_ext_params_cpy},
  gga_c_apbe_bfit_init, NULL,
  NULL, work_gga, NULL
};


static void
hyb_gga_xc_lc_wapbe_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_APBE, XC_GGA_C_APBE};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 1.0, -1.0, 0.37);

  xc_func_set_ext_params_name(p->func_aux[0], "_omega", p->hyb_omega[0]);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_wapbe = {
  XC_HYB_GGA_XC_LC_WAPBE,
  XC_EXCHANGE_CORRELATION,
  "LC hybrid based on APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  hyb_gga_xc_lc_wapbe_init, NULL,
  NULL, NULL, NULL
};

static void
hyb_gga_xc_lc_wpbeint_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_PBEINT, XC_GGA_C_PBEINT};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 1.0, -1.0, 0.52);

  xc_func_set_ext_params_name(p->func_aux[0], "_omega", p->hyb_omega[0]);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_wpbeint = {
  XC_HYB_GGA_XC_LC_WPBEINT,
  XC_EXCHANGE_CORRELATION,
  "LC hybrid based on PBEINT",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  hyb_gga_xc_lc_wpbeint_init, NULL,
  NULL, NULL, NULL
};

static void
hyb_gga_xc_lc_wsg4_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_SG4, XC_GGA_C_SG4};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 1.0, -1.0, 0.50);

  xc_func_set_ext_params_name(p->func_aux[0], "_omega", p->hyb_omega[0]);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_wsg4 = {
  XC_HYB_GGA_XC_LC_WSG4,
  XC_EXCHANGE_CORRELATION,
  "LC hybrid based on PBEINT",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  hyb_gga_xc_lc_wsg4_init, NULL,
  NULL, NULL, NULL
};
