/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_NCAP 180 /* Nearly correct asymptotic potential */

typedef struct{
  double zeta;
  double mu;
  double alpha;
  double beta;
} gga_x_ncap_params;

static void
gga_x_ncap_init(xc_func_type *p)
{
  gga_x_ncap_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ncap_params));
  params = (gga_x_ncap_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/gga_exc/gga_x_ncap.c"
#include "work_gga_new.c"

static const func_params_type ext_params[] = {
  {"_alpha", 0.34511172247159020479L, "alpha"}, /* alpha parameter */
  {"_beta", 0.018085697L, "beta"}, /* beta parameter */
  {"_mu", 0.219514973L, "mu"}, /* mu parameter */
  {"_zeta", 0.304121419L, "zeta"} /* zeta parameter */
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_ncap_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ncap_params *) (p->params);

  params->alpha = get_ext_param(p->info->ext_params, ext_params, 0);
  params->beta = get_ext_param(p->info->ext_params, ext_params, 1);
  params->mu = get_ext_param(p->info->ext_params, ext_params, 2);
  params->zeta = get_ext_param(p->info->ext_params, ext_params, 3);
}

const xc_func_info_type xc_func_info_gga_x_ncap = {
  XC_GGA_X_NCAP,
  XC_EXCHANGE,
  "Nearly correct asymptotic potential",
  XC_FAMILY_GGA,
  {&xc_ref_Carmona2019_303, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  4, ext_params, set_ext_params,
  gga_x_ncap_init, NULL,
  NULL, work_gga, NULL
};
