/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define  XC_HYB_GGA_XC_LC_BLYP 400  /* Long-range corrected BLYP */

void
xc_hyb_gga_xc_lc_blyp_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_ITYH, XC_GGA_C_LYP};
  static double funcs_coef[2];

  double gamma = 0.3;

  funcs_coef[0] = 1.0;
  funcs_coef[1] = 1.0;

  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[0], &gamma);

  xc_hyb_init_cam(p, 1.0, -1.0, gamma);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_blyp = {
  XC_HYB_GGA_XC_LC_BLYP,
  XC_EXCHANGE_CORRELATION,
  "LC version of BLYP",
  XC_FAMILY_GGA,
  {&xc_ref_Anderson2017_1656, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-9,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_lc_blyp_init, NULL,
  NULL, NULL, NULL
};
