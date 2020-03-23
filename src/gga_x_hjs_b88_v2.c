/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_HJS_B88_V2   46 /* HJS screened exchange corrected B88 version */

typedef struct{
  double omega;
} gga_x_hjs_b88_v2_params;

static void
gga_x_hjs_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hjs_b88_v2_params));
}

static const char  *omega_names[]  = {"omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.11};

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_hjs_b88_v2.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_b88_v2 = {
  XC_GGA_X_HJS_B88_V2,
  XC_EXCHANGE,
  "HJS screened exchange B88 corrected version",
  XC_FAMILY_GGA,
  {&xc_ref_Weintraub2009_754, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-6, /* densities smaller than 1e-6 yield NaNs */
  {1, omega_names, omega_desc, omega_values, set_ext_params_omega},
  gga_x_hjs_init, NULL, 
  NULL, work_gga, NULL
};
