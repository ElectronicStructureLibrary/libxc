/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_K_PC07          543 /* Perdew and Constantin 2007 */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_k_pc07.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_k_pc07 = {
  XC_MGGA_K_PC07,
  XC_KINETIC,
  "Perdew and Constantin 2007",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2007_155109, NULL, NULL, NULL, NULL},
  XC_FLAGS_DEVELOPMENT | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
