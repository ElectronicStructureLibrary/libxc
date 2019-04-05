/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_BN05   588   /* Baer and Neuhauser, gamma=1 */

#include "maple2c/lda_exc/lda_xc_bn05.c"
#include "work_lda_new.c"

const xc_func_info_type xc_func_info_lda_xc_bn05 = {
  XC_LDA_XC_BN05,
  XC_EXCHANGE_CORRELATION,
  "Baer and Neuhauser, gamma=1",
  XC_FAMILY_LDA,
  {&xc_ref_Baer2005_043002, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-24,
  0, NULL, NULL,
  NULL, NULL,
  work_lda, NULL, NULL
};
