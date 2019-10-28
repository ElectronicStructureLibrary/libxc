/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_OL2          183 /* Exchange form based on Ou-Yang and Levy v.2 */

typedef struct{
  double aa, bb, cc;
} gga_x_ol2_params;

static void 
gga_x_ol2_init(xc_func_type *p)
{
  gga_x_ol2_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ol2_params));
  params = (gga_x_ol2_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/gga_exc/gga_x_ol2.c"
#include "work_gga.c"

static const func_params_type ext_params[] = {
  {"_aa", 0.09564574034649151285038696952714226444963L, "aa"}, /* M_CBRT2*0.07064/X_FACTOR_C */
  {"_bb", 0.09564574034649151285038696952714226444963L, "bb"}, /* M_CBRT2*0.07064/X_FACTOR_C */
  {"_cc", 4.098833606342553442039881031486386917472L,   "cc"}  /* M_CBRT2*M_CBRT2*0.07064*34.0135/X_FACTOR_C */
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_ol2_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ol2_params *) (p->params);

  params->aa = get_ext_param(p->info->ext_params, ext_params, 0);
  params->bb = get_ext_param(p->info->ext_params, ext_params, 1);
  params->cc = get_ext_param(p->info->ext_params, ext_params, 2);
}

const xc_func_info_type xc_func_info_gga_x_ol2 = {
  XC_GGA_X_OL2,
  XC_EXCHANGE,
  "Exchange form based on Ou-Yang and Levy v.2",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, &xc_ref_OuYang1991_379, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-26,
  3, ext_params, set_ext_params,
  gga_x_ol2_init, NULL, 
  NULL, work_gga, NULL
};
