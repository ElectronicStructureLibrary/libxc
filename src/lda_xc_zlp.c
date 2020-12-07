/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_XC_ZLP     43   /* Zhao, Levy & Parr, Eq. (20)  */

typedef struct {
  double zlp_a0, zlp_k;
} lda_xc_zlp_params;

#define N_PAR 2
static const char *names[N_PAR] = {"_a0", "_k"};
static const char *desc[N_PAR] = {"a0 parameter", "k parameter"};
static const double par_zlp[N_PAR] = {0.93222, 9.47362e-3};

static void
lda_xc_zlp_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_xc_zlp_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_xc_zlp.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_zlp = {
  XC_LDA_XC_ZLP,
  XC_EXCHANGE_CORRELATION,
  "Zhao, Levy & Parr, Eq. (20)",
  XC_FAMILY_LDA,
  {&xc_ref_Zhao1993_918, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  {N_PAR, names, desc, par_zlp, set_ext_params_cpy},
  lda_xc_zlp_init, NULL,
  work_lda, NULL, NULL
};
