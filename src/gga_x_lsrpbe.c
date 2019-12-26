/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2018      Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_LSRPBE  169 /* PW91-like modification of RPBE */

typedef struct{
  double kappa;
  double mu;
  double alpha;
} gga_x_lsrpbe_params;

static void 
gga_x_lsrpbe_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_lsrpbe_params));
}

static const func_params_type ext_params[] = {
  {"_kappa", 0.8040, "Asymptotic value of the enhancement function"},
  {"_mu",    MU_PBE, "Coefficient of the 2nd order expansion of the full Lsrpbe functional"},
  {"_alpha", 0.00680892, "Exponent that should satisfy the PW91 criterion"}
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_lsrpbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_lsrpbe_params *) (p->params);

  params->kappa = get_ext_param(p->info->ext_params, ext_params, 0);
  params->mu    = get_ext_param(p->info->ext_params, ext_params, 1);
  params->alpha = get_ext_param(p->info->ext_params, ext_params, 2);

  /* adapt used mu value to yield wanted mu near origin (eq 9) */
  params-> mu += params->alpha*(1.0 + params->kappa);

  printf("mu = %.10f\n",params->mu);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_lsrpbe.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_lsrpbe = {
  XC_GGA_X_LSRPBE,
  XC_EXCHANGE,
  "lsRPBE, a PW91-like modification of RPBE",
  XC_FAMILY_GGA,
  {&xc_ref_PachecoKato2016_268, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  3, ext_params, set_ext_params,
  gga_x_lsrpbe_init, NULL, 
  NULL, work_gga, NULL
};
