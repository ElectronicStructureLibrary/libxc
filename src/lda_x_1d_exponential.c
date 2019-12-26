/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_X_1D_EXPONENTIAL  600 /* Exchange in 1D for an exponentially screened interaction */

typedef struct{
  double beta;         /* screening parameter beta */
} lda_x_1d_exponential_params;

static void 
lda_x_1d_exponential_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_x_1d_exponential_params));
}

static inline double FT_inter(double x)
{
  double x2 = x*x;
  return expint_e1(x2)*exp(x2);
}

static void func1(double *x, int n, void *dummy)
{
  int ii;
  
  for(ii=0; ii<n; ii++)
    x[ii] = FT_inter(x[ii]);
}

static void func2(double *x, int n, void *dummy)
{
  int ii;
  
  for(ii=0; ii<n; ii++)
    x[ii] = x[ii]*FT_inter(x[ii]);
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_x_1d_exponential.c"
#include "work_lda.c"

static const func_params_type ext_params[] = {
  {"beta", 1.0, "Screening parameter"}
};
static void 

set_ext_params(xc_func_type *p, const double *ext_params)
{
  lda_x_1d_exponential_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_x_1d_exponential_params *)(p->params);

  params->beta = get_ext_param(p->info->ext_params, ext_params, 0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_x_1d_exponential = {
  XC_LDA_X_1D_EXPONENTIAL,
  XC_EXCHANGE,
  "Exchange in 1D for an exponentially screened interaction",
  XC_FAMILY_LDA,
  {&xc_ref_Helbig2011_032503, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | MAPLE2C_FLAGS,
  1e-26,
  1, ext_params, set_ext_params,
  lda_x_1d_exponential_init, NULL,
  work_lda, NULL, NULL
};
