/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_RTPSS          299 /* Revised TPSS exchange by Garza, Bell and Head-Gordon */

typedef struct {
  double b, c, e, kappa, mu;
} mgga_x_rtpss_params;

static void 
mgga_x_rtpss_init(xc_func_type *p)
{
  mgga_x_rtpss_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_rtpss_params));
  params = (mgga_x_rtpss_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_RTPSS:
    /* default set by set_ext_params */
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_rtpss\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_b",      0.40,    "b"},
  {"_c",      1.59096, "c"},
  {"_e",      1.537,   "e"},
  {"_kappa",  0.8040,  "Asymptotic value of the enhancement function"},
  {"_mu",     0.21951, "Coefficient of the 2nd order expansion"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_rtpss_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_rtpss_params *) (p->params);

  params->b      = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c      = get_ext_param(p->info->ext_params, ext_params, 1);
  params->e      = get_ext_param(p->info->ext_params, ext_params, 2);
  params->kappa  = get_ext_param(p->info->ext_params, ext_params, 3);
  params->mu     = get_ext_param(p->info->ext_params, ext_params, 4);
}

#include "maple2c/mgga_x_rtpss.c"

#define func xc_mgga_x_rtpss_enhance
#include "work_mgga_x.c"

const xc_func_info_type xc_func_info_mgga_x_rtpss = {
  XC_MGGA_X_RTPSS,
  XC_EXCHANGE,
  "TPSS for surface adsorption",
  XC_FAMILY_MGGA,
  {&xc_ref_Garza2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23,
  5, ext_params, set_ext_params,
  mgga_x_rtpss_init, NULL, 
  NULL, NULL, work_mgga_x,
};
