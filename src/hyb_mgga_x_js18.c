/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_MGGA_X_JS18       705 /* a screened version of TM */

static void
hyb_mgga_x_js18_init(xc_func_type *p)
{
  xc_hyb_init_sr(p, 0.33, 0.1);
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_x_js18.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_js18 = {
  XC_HYB_MGGA_X_JS18,
  XC_EXCHANGE,
  "JS18",
  XC_FAMILY_MGGA,
  {&xc_ref_Jana2018_8999, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  hyb_mgga_x_js18_init, NULL,
  NULL, NULL, work_mgga
};
