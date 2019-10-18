/*
 Copyright (C) 2017 Miguel Marques, Mario Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_GX          575 /* GX functional of Loos */

typedef struct{
  double c0, c1, alphainf;
} mgga_x_gx_params;

static void
mgga_x_gx_init(xc_func_type *p)
{
  mgga_x_gx_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_gx_params));
  params = (mgga_x_gx_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/mgga_exc/mgga_x_gx.c"
#include "work_mgga.c"

static const func_params_type ext_params[] = {
  {"_c0", 0.827411L, "c0"}, /* c_0 */
  {"_c1", -0.643560L, "c1"}, /* c_1 */
  {"_alphainf", 0.852, "alphainf"}, /* \alpha_\infty */
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_gx_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_gx_params *) (p->params);

  params->c0 = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c1 = get_ext_param(p->info->ext_params, ext_params, 1);
  params->alphainf = get_ext_param(p->info->ext_params, ext_params, 2);
}

const xc_func_info_type xc_func_info_mgga_x_gx = {
  XC_MGGA_X_GX,
  XC_EXCHANGE,
  "GX functional of Loos",
  XC_FAMILY_MGGA,
  {&xc_ref_Loos2017_114108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-20,
  3, ext_params, set_ext_params,
  mgga_x_gx_init, NULL,
  NULL, NULL, work_mgga,
};
