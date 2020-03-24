/*
 Copyright (C) 2013 Rolf Wuerdemann, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_SFAT 530 /* short-range recipe for PBE functional */
/* see
   Savin, Flad, Int. J. Quant. Chem. Vol. 56, 327-332 (1995)
   Akinaga, Ten-no, Chem. Phys. Lett. 462 (2008) 348-351
*/

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_sfat.c"
#include "work_gga.c"

static const char  *omega_names[]  = {"omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.44};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_sfat = {
  XC_GGA_X_SFAT,
  XC_EXCHANGE,
  "Short-range recipe for B88 functional - Yukawa",
  XC_FAMILY_GGA,
  {&xc_ref_Savin1995_327, &xc_ref_Akinaga2008_348, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-19,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  NULL, NULL,
  NULL, work_gga, NULL
};

