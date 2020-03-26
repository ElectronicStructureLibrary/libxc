/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_ITYH 529 /* short-range recipe B88 functionals - erf */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ityh.c"
#include "work_gga.c"

static const char  *omega_names[]  = {"omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.2};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ityh = {
  XC_GGA_X_ITYH,
  XC_EXCHANGE,
  "Short-range recipe for B88 functional - erf",
  XC_FAMILY_GGA,
  {&xc_ref_Iikura2001_3540, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-8,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  NULL, NULL, 
  NULL, work_gga, NULL
};
