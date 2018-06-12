/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_ERF   546   /* Attenuated exchange LDA (erf) */

static const func_params_type ext_params[] = {
  {"omega",  0.3, "screening parameter"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);

  p->cam_omega = get_ext_param(p->info->ext_params, ext_params, 0);
}


#include "maple2c/lda_exc/lda_x_erf.c"
#include "work_lda_new.c"

const xc_func_info_type xc_func_info_lda_x_erf = {
  XC_LDA_X_ERF,
  XC_EXCHANGE,
  "Attenuated exchange LDA (erf)",
  XC_FAMILY_LDA,
  {&xc_ref_Gill1996_1005, &xc_ref_Toulouse2004_1047, &xc_ref_Tawada2004_8425, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-13,
  1, ext_params, set_ext_params,
  NULL, NULL, 
  work_lda, NULL, NULL
};
