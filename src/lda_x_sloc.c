/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_SLOC   692   /* simple local model for Slater potential */

#include "maple2c/lda_exc/lda_x_sloc.c"
#include "work_lda_new.c"

const xc_func_info_type xc_func_info_lda_x_sloc = {
  XC_LDA_X_SLOC,
  XC_EXCHANGE,
  "simple local model for Slater potential",
  XC_FAMILY_LDA,
  {&xc_ref_Finzel2017_40, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  0, NULL, NULL,
  NULL, NULL, 
  work_lda, NULL, NULL
};
