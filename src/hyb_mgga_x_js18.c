/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_JS18       705 /* a screened version of TM */

static void
hyb_mgga_x_js18_init(xc_func_type *p)
{
  p->cam_omega =  0.33;
  p->cam_alpha =  0.1;
}

const xc_func_info_type xc_func_info_hyb_gga_xc_js18 = {
  XC_HYB_GGA_XC_JS18,
  XC_EXCHANGE,
  "JS18",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Heyd2003_8207, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_x_js18_init, NULL,
  NULL, NULL, NULL
};
