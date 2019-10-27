/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_K_GDS08_WORKER 100001   /* Combined analytical theory with Monte Carlo sampling */

typedef struct {
  double A, B, C;
} lda_k_gds08_params;

static void 
lda_k_gds08_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_k_gds08_params));
}

static func_params_type ext_params[] = {
  {"_A", 0.860, "linear term"},
  {"_B", 0.224, "term proportional to the logarithm of the density"},
  {"_C", 0.0,   "term proportional to the square of the logarithm"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_k_gds08_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_k_gds08_params *) (p->params);

  params->A = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B = get_ext_param(p->info->ext_params, ext_params, 1);
  params->C = get_ext_param(p->info->ext_params, ext_params, 2);
}

#include "maple2c/lda_exc/lda_k_gds08_worker.c"
#include "work_lda.c"

const xc_func_info_type xc_func_info_lda_k_gds08_worker = {
  XC_LDA_K_GDS08_WORKER,
  XC_KINETIC,
  "Combined analytical theory with Monte Carlo sampling",
  XC_FAMILY_LDA,
  {&xc_ref_Ghiringhelli2008_073104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  3, ext_params, set_ext_params,
  lda_k_gds08_init, NULL,
  work_lda, NULL, NULL
};
