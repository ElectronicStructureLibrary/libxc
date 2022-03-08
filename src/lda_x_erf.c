/*
 Copyright (C) 2017 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_ERF         546 /* Short-range LDA exchange with error function kernel (erfc) */
#define  XC_HYB_LDA_X_ERF    744 /* Hybrid Functional with short-ranged LDA exchange (erfc)*/

#include "maple2c/lda_exc/lda_x_erf.c"
#include "work_lda.c"

static void
xc_lda_x_erf_init(xc_func_type *p)
{
  assert(p!=NULL);
  xc_hyb_init_hybrid(p, 0.0);
  p->hyb_type[0] = XC_HYB_NONE;
}

static const char  *omega_names[]  = {"_omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.3};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_erf = {
  XC_LDA_X_ERF,
  XC_EXCHANGE,
  "Short-range LDA exchange with error function kernel (erfc)",
  XC_FAMILY_LDA,
  {&xc_ref_Gill1996_1005, &xc_ref_Toulouse2004_1047, &xc_ref_Tawada2004_8425, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_omega},
  xc_lda_x_erf_init, NULL,
  &work_lda, NULL, NULL
};

static void
xc_hyb_lda_x_erf_init(xc_func_type *p)
{
  assert(p!=NULL);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

const xc_func_info_type xc_func_info_hyb_lda_x_erf = {
  XC_HYB_LDA_X_ERF,
  XC_EXCHANGE,
  "Hybrid Exchange functional with short-range LDA exchange (erf)",
  XC_FAMILY_LDA,
  {&xc_ref_Gill1996_1005, &xc_ref_Toulouse2004_1047, &xc_ref_Tawada2004_8425, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {1, omega_names, omega_desc, omega_values, set_ext_params_cpy_lc},
  xc_hyb_lda_x_erf_init, NULL,
  work_lda, NULL, NULL
};

