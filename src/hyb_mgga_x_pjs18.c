/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_MGGA_X_PJS18       706 /* a screened version of TM */

static void
hyb_mgga_x_pjs18_init(xc_func_type *p)
{
  xc_hyb_init_cam(p, 0.33, 1.0, -1.0);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_x_pjs18.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_pjs18 = {
  XC_HYB_MGGA_X_PJS18,
  XC_EXCHANGE,
  "PJS18",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Patra2018_8991, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  hyb_mgga_x_pjs18_init, NULL,
  NULL, NULL, work_mgga
};
