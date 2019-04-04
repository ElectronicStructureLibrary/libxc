/*
 Copyright (C) 2018 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_MVSB          302 /* MVSBeta exchange of Furness and Sun */
#define XC_MGGA_X_MVSBS         303 /* MVSBeta* exchange of Furness and Sun */

typedef struct {
  double e1, c1, k0, b;
} mgga_x_mvsb_params;

static void
mgga_x_mvsb_init(xc_func_type *p)
{
  mgga_x_mvsb_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_mvsb_params));
  params = (mgga_x_mvsb_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_MVSB:
    /* set by set_ext_params */
    break;
  case XC_MGGA_X_MVSBS:
    params->e1 = -2.3800;
    params->c1 =  6.3783;
    params->k0 =  0.174;
    params->b  =  0.0233;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_mvs\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_e1", -1.6665, "e1 parameter"},
  {"_c1", 7.8393, "c1 parameter"},
  {"_k0", 0.174, "k0 parameter"},
  {"_b", 0.0233, "b parameter"}
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_mvsb_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_mvsb_params *) (p->params);

  params->e1 = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c1 = get_ext_param(p->info->ext_params, ext_params, 1);
  params->k0 = get_ext_param(p->info->ext_params, ext_params, 2);
  params->b  = get_ext_param(p->info->ext_params, ext_params, 3);
}


#include "maple2c/mgga_exc/mgga_x_mvsb.c"
#include "work_mgga_new.c"

const xc_func_info_type xc_func_info_mgga_x_mvsb = {
  XC_MGGA_X_MVSB,
  XC_EXCHANGE,
  "MVSbeta exchange by Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_ALL,
  1e-23,
  4, ext_params, set_ext_params,
  mgga_x_mvsb_init, NULL,
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_mgga_x_mvsbs = {
  XC_MGGA_X_MVSBS,
  XC_EXCHANGE,
  "MVSbeta* exchange by Furness and Sun",
  XC_FAMILY_MGGA,
  {&xc_ref_Furness2018, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_ALL,
  1e-23,
  0, NULL, NULL,
  mgga_x_mvsb_init, NULL,
  NULL, NULL, work_mgga,
};
