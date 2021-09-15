/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_APBE_BFIT    720 /* APBE with a fit of beta */
#define XC_GGA_C_PBEINT_BFIT  721 /* PBEINT with a fit of beta */

typedef struct{
  double mux, omega;
} gga_c_apbe_bfit_params;


static void gga_c_apbe_bfit_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_apbe_bfit_params));
}

#define APBE_BFIT_N_PAR 2
static const char  *apbe_bfit_names[APBE_BFIT_N_PAR]  = {"_mux", "_omega"};
static const char  *apbe_bfit_desc[APBE_BFIT_N_PAR]   = {
  "mu of the corresponding exchange functional",
  "range separation parameter"};
static const double apbe_bfit_values[APBE_BFIT_N_PAR] = 
  {0.26, 0.37};
static const double pbeint_bfit_values[APBE_BFIT_N_PAR] = 
  {10./81., 0.52};

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
