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
  p->cam_omega =  0.33;
  p->cam_beta  =  0.1;
}

#include "maple2c/mgga_exc/hyb_mgga_x_pjs18.c"
#include "work_mgga_new.c"

const xc_func_info_type xc_func_info_hyb_mgga_x_pjs18 = {
  XC_HYB_MGGA_X_PJS18,
  XC_EXCHANGE,
  "PJS18",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Patra2018_8991, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_x_pjs18_init, NULL,
  NULL, NULL, work_mgga
};
