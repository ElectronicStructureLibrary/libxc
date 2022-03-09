/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HGGA_X_JSE_B88          735  /* Local hyrbid functional of Jaramillo et al */

#include "maple2c/hgga_exc/hgga_x_jse_b88.c"
#include "work_hgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hgga_x_jse_b88 = {
  XC_HGGA_X_JSE_B88,
  XC_EXCHANGE,
  "Local hyrbid functional of Jaramillo et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Jaramillo2003_1068, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, NULL, NULL, &work_hgga,
};

