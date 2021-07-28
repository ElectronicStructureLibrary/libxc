/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_LDA_XC_BN05   588   /* Baer and Neuhauser, gamma=1 */

typedef struct{
  double omega;
} hyb_lda_xc_bn05_params;

#include "decl_lda.h"
#include "maple2c/lda_exc/hyb_lda_xc_bn05.c"
#include "work_lda.c"

#define N_PAR 1
static const char  *names[N_PAR]  = {"_omega"};
static const char  *desc[N_PAR]   = {
  "Range separation parameter"
};

static const double par_bn05[N_PAR] = {1.0};

static void
hyb_lda_xc_bn05_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(hyb_lda_xc_bn05_params));

  xc_hyb_init_camy(p, 1.0, -1.0, par_bn05[0]);
}

static void
bn05_set_ext_params(xc_func_type *p, const double *ext_params)
{
  hyb_lda_xc_bn05_params *params;

  assert(p->params != NULL);
  params = (hyb_lda_xc_bn05_params * )(p->params);

  set_ext_params_cpy(p, ext_params);
  p->hyb_params[1].sr.omega = params->omega;
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_lda_xc_bn05 = {
  XC_HYB_LDA_XC_BN05,
  XC_EXCHANGE_CORRELATION,
  "Baer and Neuhauser, gamma=1",
  XC_FAMILY_LDA,
  {&xc_ref_Baer2005_043002, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR, names, desc, par_bn05, bn05_set_ext_params},
  hyb_lda_xc_bn05_init, NULL,
  work_lda, NULL, NULL
};
