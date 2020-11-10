/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_R2SCAN         498 /* Re-regularized SCAN correlation */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_r2scan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_r2scan = {
  XC_MGGA_C_R2SCAN,
  XC_CORRELATION,
  "Re-regularized SCAN correlation by Furness et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2020_8208, &xc_ref_Furness2020_9248, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga,
};
