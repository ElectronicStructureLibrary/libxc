/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_BR89_EXPLICIT    586 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8 */
#define XC_MGGA_X_BR89_EXPLICIT_1  602 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0 */

typedef struct{
  double gamma;
} mgga_x_br89_params;

static const mgga_x_br89_params par_one = {1.0};

static void
mgga_x_br89_init(xc_func_type *p)
{
  mgga_x_br89_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_br89_params));
  params = (mgga_x_br89_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_BR89_EXPLICIT:
    /* default set by set_ext_params */
    break;
  case XC_MGGA_X_BR89_EXPLICIT_1:
    memcpy(params, &par_one, sizeof(mgga_x_br89_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_br89_explicit\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_gamma", 0.8, "gamma"},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_br89_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_br89_params *) (p->params);

  params->gamma = get_ext_param(p->info->ext_params, ext_params, 0);
}

#include "maple2c/mgga_x_br89_explicit.c"

#define func maple2c_func
#include "work_mgga_x.c"

const xc_func_info_type xc_func_info_mgga_x_br89_explicit = {
  XC_MGGA_X_BR89_EXPLICIT,
  XC_EXCHANGE,
  "Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1989_3761, &xc_ref_Proynov2008_103, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1.0e-12,
  1, ext_params, set_ext_params,
  mgga_x_br89_init, NULL,
  NULL, NULL, work_mgga_x
};

const xc_func_info_type xc_func_info_mgga_x_br89_explicit_1 = {
  XC_MGGA_X_BR89_EXPLICIT_1,
  XC_EXCHANGE,
  "Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1989_3761, &xc_ref_Proynov2008_103, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_ALL,
  1.0e-12,
  0, NULL, NULL,
  mgga_x_br89_init, NULL,
  NULL, NULL, work_mgga_x
};
