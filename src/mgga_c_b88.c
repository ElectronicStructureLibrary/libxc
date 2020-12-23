/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_B88          571 /* Meta-GGA correlation by Becke */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_b88.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_b88 = {
  XC_MGGA_C_B88,
  XC_CORRELATION,
  "Meta-GGA correlation by Becke",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1988_1053, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-14,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
