/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_RSCANL         706 /* Regularized SCAN correlation */

#include "mgga_k_pcopt.c"
#include "maple2c/mgga_exc/mgga_c_rscan.c"
#include "work_mgga_deorb.c"

const xc_func_info_type xc_func_info_mgga_c_rscanl = {
  XC_MGGA_C_RSCANL,
  XC_CORRELATION,
  "Deorbitalized and regularized SCAN (rSCAN-L) correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Bartok2019_161101, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_HAVE_EXC,
  1e-20,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, NULL, work_mgga_deorb,
};
