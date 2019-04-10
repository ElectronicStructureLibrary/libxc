/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_TIH   599   /* Neural network LDA from Tozer et al */

#include "maple2c/lda_vxc/lda_xc_tih.c"
#define XC_NO_EXC
#include "work_lda_new.c"

const xc_func_info_type xc_func_info_lda_xc_tih = {
  XC_LDA_XC_TIH,
  XC_EXCHANGE_CORRELATION,
  "Neural network LDA from Tozer et al",
  XC_FAMILY_LDA,
  {&xc_ref_Tozer1996_9200, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC,
  5e-24,
  0, NULL, NULL,
  NULL, NULL,
  work_lda, NULL, NULL
};
