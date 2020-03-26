/*
 Copyright (C) 2019 M.A.L. Marques
                    Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_TASK      707 /* TASK exchange of Aschebrock and Kuemmel */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_task.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_task = {
  XC_MGGA_X_TASK,
  XC_EXCHANGE,
  "TASK exchange of Aschebrock and Kuemmel",
  XC_FAMILY_MGGA,
  {&xc_ref_Aschebrock2019_033082, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, work_mgga
};
