/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* these functionals is for the soft-Coulomb interaction with alpha=1 */

#define XC_LDA_XC_1D_EHWLRG_1     536 /* LDA constructed from slab-like systems of 1 electron  */
#define XC_LDA_XC_1D_EHWLRG_2     537 /* LDA constructed from slab-like systems of 2 electrons */
#define XC_LDA_XC_1D_EHWLRG_3     538 /* LDA constructed from slab-like systems of 3 electrons */

typedef struct {
  double ehwlrg_alpha;
  double ehwlrg_a1, ehwlrg_a2, ehwlrg_a3;
} lda_xc_1d_ehwlrg_params;

#define N_PAR 4
static const char *names[N_PAR] = {"_alpha", "_a1", "_a2", "_a3"};
static const char *desc[N_PAR] = {"alpha", "a1", "a2", "a3"};

static const double par_ehwlrg_1[N_PAR] = {0.638, -0.803, 0.82, -0.47};
static const double par_ehwlrg_2[N_PAR] = {0.604, -0.74, 0.68, -0.38};
static const double par_ehwlrg_3[N_PAR] = {0.61, -0.77, 0.79, -0.48};
  
static void
lda_xc_1d_ehwlrg_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_xc_1d_ehwlrg_params));
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_xc_1d_ehwlrg.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_1d_ehwlrg_1 = {
  XC_LDA_XC_1D_EHWLRG_1,
  XC_EXCHANGE_CORRELATION,
  "LDA constructed from slab-like systems of 1 electron",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-32,
  {N_PAR, names, desc, par_ehwlrg_1, set_ext_params_cpy},
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_1d_ehwlrg_2 = {
  XC_LDA_XC_1D_EHWLRG_2,
  XC_EXCHANGE_CORRELATION,
  "LDA constructed from slab-like systems of 2 electrons",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-32,
  {N_PAR, names, desc, par_ehwlrg_2, set_ext_params_cpy},
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_1d_ehwlrg_3 = {
  XC_LDA_XC_1D_EHWLRG_3,
  XC_EXCHANGE_CORRELATION,
  "LDA constructed from slab-like systems of 3 electrons",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-32,
  {N_PAR, names, desc, par_ehwlrg_3, set_ext_params_cpy},
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};
