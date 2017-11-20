/*
 Copyright (C) 2006-2016 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_C_CHACHIYO  287   /* Chachiyo simple 2 parameter correlation   */

#include "maple2c/lda_c_chachiyo.c"

#define func maple2c_func
#include "work_lda.c"

const xc_func_info_type xc_func_info_lda_c_chachiyo = {
  XC_LDA_C_CHACHIYO,
  XC_CORRELATION,
  "Chachiyo simple 2 parameter correlation",
  XC_FAMILY_LDA,
  {&xc_ref_Chachiyo2016_021101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
  NULL,
  NULL
};
