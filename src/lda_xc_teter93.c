/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_XC_TETER93     20   /* Teter 93 parametrization                */

typedef struct {
  double teter_a[4], teter_ap[4];
  double teter_b[4], teter_bp[4];
} lda_xc_teter_params;

#define N_PAR 16
static const char *names[N_PAR] = {
  "_a0", "_a1", "_a2", "_a3", "_ap0", "_ap1", "_ap2", "_ap3",
  "_b0", "_b1", "_b2", "_b3", "_bp0", "_bp1", "_bp2", "_bp3"
};
static const char *desc[N_PAR] = {
  "a0", "a1", "a2", "a3", "ap0", "ap1", "ap2", "ap3",
  "b0", "b1", "b2", "b3", "bp0", "bp1", "bp2", "bp3"
};

static const double par_teter[N_PAR] = {
  0.4581652932831429, 2.217058676663745,  0.7405551735357053, 0.01968227878617998,
  0.119086804055547,  0.6157402568883345, 0.1574201515892867, 0.003532336663397157,
  1.0000000000000000, 4.504130959426697,  1.110667363742916,  0.02359291751427506,
  0.000000000000000,  0.2673612973836267, 0.2052004607777787, 0.004200005045691381
};


static void
lda_xc_teter_init(xc_func_type *p) {
  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_xc_teter_params));
}


#include "decl_lda.h"
#include "maple2c/lda_exc/lda_xc_teter93.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_xc_teter93 = {
  XC_LDA_XC_TETER93,
  XC_EXCHANGE_CORRELATION,
  "Teter 93",
  XC_FAMILY_LDA,
  {&xc_ref_Goedecker1996_1703, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_teter, set_ext_params_cpy},
  lda_xc_teter_init, NULL,
  work_lda, NULL, NULL
};
