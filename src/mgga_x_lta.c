/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_LTA          201 /* Local tau approximation of Ernzerhof & Scuseria */
#define XC_MGGA_X_TLDA         685 /* LDA-type exchange with tau-dependent potential */

typedef struct{
  double power;
} mgga_x_lta_params;

static void 
mgga_x_lta_init(xc_func_type *p)
{
  mgga_x_lta_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_lta_params));
  params = (mgga_x_lta_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_LTA:
    /* default values set by set_ext_params */
    break;
  case XC_MGGA_X_TLDA:
    params->power = 1.0/5.0;
    break;
  default:     
    fprintf(stderr, "Internal error in mgga_x_lta\n");
    exit(1);
  } 
}

static func_params_type ext_params[] = {
  {"_power", 4.0/5.0, "power of t"},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_lta_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_lta_params *) (p->params);

  params->power  = get_ext_param(p->info->ext_params, ext_params, 0);
}

#include "maple2c/mgga_exc/mgga_x_lta.c"
#include "work_mgga.c"

const xc_func_info_type xc_func_info_mgga_x_lta = {
  XC_MGGA_X_LTA,
  XC_EXCHANGE,
  "Local tau approximation",
  XC_FAMILY_MGGA,
  {&xc_ref_Ernzerhof1999_911, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-23,
  1, ext_params, set_ext_params,
  mgga_x_lta_init, NULL,
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_mgga_x_tlda = {
  XC_MGGA_X_TLDA,
  XC_EXCHANGE,
  "LDA-type exchange with tau-dependent potential",
  XC_FAMILY_MGGA,
  {&xc_ref_Eich2014_224107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-23,
  0, NULL, NULL,
  mgga_x_lta_init, NULL,
  NULL, NULL, work_mgga,
};
