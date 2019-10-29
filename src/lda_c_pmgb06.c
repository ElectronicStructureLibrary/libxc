/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_PMGB06   590   /* Long-range LDA correlation functional */

static const func_params_type ext_params[] = {
  {"omega",  0.3, "screening parameter"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);

  p->cam_omega = get_ext_param(p->info->ext_params, ext_params, 0);
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_pmgb06.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_pmgb06 = {
  XC_LDA_C_PMGB06,
  XC_EXCHANGE,
  "Long-range LDA correlation functional",
  XC_FAMILY_LDA,
  {&xc_ref_Paziani2006_155111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS | XC_FLAGS_DEVELOPMENT,
  1e-13,
  1, ext_params, set_ext_params,
  NULL, NULL, 
  work_lda, NULL, NULL
};
