/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_BC95          240 /* Becke correlation 95 */

typedef struct{
  double css, copp;
} mgga_c_bc95_params;


static void 
mgga_c_bc95_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_bc95_params));
}

static const func_params_type ext_params[] = {
  {"_css",  0.038,  "Parallel spin"},
  {"_copp", 0.0031, "Opposite spin"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_c_bc95_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_c_bc95_params *) (p->params);

  params->css  = get_ext_param(p->info->ext_params, ext_params, 0);
  params->copp = get_ext_param(p->info->ext_params, ext_params, 1);
}

#include "maple2c/mgga_exc/mgga_c_bc95.c"
#include "work_mgga.c"

const xc_func_info_type xc_func_info_mgga_c_bc95 = {
  XC_MGGA_C_BC95,
  XC_CORRELATION,
  "Becke correlation 95",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1996_1040, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-23,
  2, ext_params, set_ext_params,
  mgga_c_bc95_init, NULL,
  NULL, NULL, work_mgga,
};

