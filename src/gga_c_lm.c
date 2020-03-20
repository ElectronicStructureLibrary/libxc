/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/************************************************************************
  This functional is provided for historical reasons.
  It was one of the first GGAs that ever appeared.
************************************************************************/

#define XC_GGA_C_LM          137 /* Langreth and Mehl correlation          */

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_lm.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_lm = {
  XC_GGA_C_LM,
  XC_CORRELATION,
  "Langreth & Mehl",
  XC_FAMILY_GGA,
  {&xc_ref_Langreth1981_446, &xc_ref_Hu1985_391, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-12,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL,
  NULL, work_gga, NULL
};

