/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_XC_BN05   588   /* Baer and Neuhauser, gamma=1 */

typedef struct {
  double bn05_A, bn05_C0, bn05_C1;
} lda_xc_bn05_params;

#define N_PAR 3
static const char *names[N_PAR] = {"_A", "_C0", "_C1"};
static const char *desc[N_PAR] = {"A", "C0", "C1"};
static const double par_bn05[N_PAR] = {3.4602, 3.2, -0.9};

static void
lda_xc_bn05_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_xc_bn05_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_xc_bn05.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_bn05 = {
  XC_LDA_XC_BN05,
  XC_EXCHANGE_CORRELATION,
  "Baer and Neuhauser, gamma=1",
  XC_FAMILY_LDA,
  {&xc_ref_Baer2005_043002, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_bn05, set_ext_params_cpy},
  lda_xc_bn05_init, NULL,
  work_lda, NULL, NULL
};
