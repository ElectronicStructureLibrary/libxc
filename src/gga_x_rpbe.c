/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_RPBE  117 /* Hammer, Hansen & Norskov (PBE-like) */


typedef struct{
  double rpbe_kappa, rpbe_mu;
} gga_x_rpbe_params;


static void 
gga_x_rpbe_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_rpbe_params));
}

static const func_params_type ext_params[] = {
  {"_kappa", 0.8040, "Asymptotic value of the enhancement function"},
  {"_mu",    MU_PBE, "Coefficient of the 2nd order expansion"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_rpbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_rpbe_params *) (p->params);

  params->rpbe_kappa = get_ext_param(p->info->ext_params, ext_params, 0);
  params->rpbe_mu    = get_ext_param(p->info->ext_params, ext_params, 1);
}

#include "maple2c/gga_exc/gga_x_rpbe.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_rpbe = {
  XC_GGA_X_RPBE,
  XC_EXCHANGE,
  "Hammer, Hansen, and Norskov",
  XC_FAMILY_GGA,
  {&xc_ref_Hammer1999_7413, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  2, ext_params, set_ext_params,
  gga_x_rpbe_init, NULL, 
  NULL, work_gga, NULL
};
