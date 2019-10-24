/*
 Copyright (C) 2017 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_X_SLOC   692   /* simple local model for Slater potential */

typedef struct{
  double a;       /* prefactor */
  double b;       /* exponent */
} lda_x_sloc_params;

static void
lda_x_sloc_init(xc_func_type *p)
{
  lda_x_sloc_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_x_sloc_params));
  params = (lda_x_sloc_params *) (p->params);

  /* default set by set_ext_params */
}

#include "maple2c/lda_exc/lda_x_sloc.c"
#include "work_lda.c"

static const func_params_type ext_params_sloc[] = {
  {"_a", 1.67, "Prefactor"},
  {"_b",  0.3, "Exponent"},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_x_sloc_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_x_sloc_params *) (p->params);

  params->a = get_ext_param(p->info->ext_params, ext_params, 0);
  params->b = get_ext_param(p->info->ext_params, ext_params, 1);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_sloc = {
  XC_LDA_X_SLOC,
  XC_EXCHANGE,
  "simple local model for Slater potential",
  XC_FAMILY_LDA,
  {&xc_ref_Finzel2017_40, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  2, ext_params_sloc, set_ext_params,
  lda_x_sloc_init, NULL,
  work_lda, NULL, NULL
};
