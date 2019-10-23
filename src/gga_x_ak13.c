/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019      Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_AK13  56 /* Armiento & Kuemmel 2013 */

typedef struct{
  double B1, B2;
} gga_x_ak13_params;

static void
gga_x_ak13_init(xc_func_type *p)
{
  gga_x_ak13_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ak13_params));
  params = (gga_x_ak13_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/gga_exc/gga_x_ak13.c"
#include "work_gga.c"

static const func_params_type ext_params[] = {
  {"_B1",  1.74959015598863046792081721182, "B1"}, /* 3*muGE/5 + 8 pi/15 */
  {"_B2", -1.62613336586517367779736042170, "B2"} /* muGE - B1 */
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_ak13_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ak13_params *) (p->params);

  params->B1 = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B2 = get_ext_param(p->info->ext_params, ext_params, 1);
}

double xc_gga_ak13_get_asymptotic (double homo)
{
  const double params[2] = {1.74959015598863046792081721182, -1.62613336586517367779736042170};
  return xc_gga_ak13_pars_get_asymptotic(homo, params);
}

double xc_gga_ak13_pars_get_asymptotic (double homo, const double *ext_params)
{
  double Qx, aa, aa2, factor;
  double ak13_B1, ak13_B2;

  ak13_B1 = ext_params[0];
  ak13_B2 = ext_params[1];

  Qx = sqrt(2.0)*ak13_B1/(3.0*CBRT(3.0*M_PI*M_PI));

  aa  = X_FACTOR_C*Qx;
  aa2 = aa*aa;

  factor = (homo < 0.0) ? -1.0 : 1.0;

  return (aa2/2.0)*(1.0 + factor*sqrt(1.0 - 4.0*homo/aa2));
}

const xc_func_info_type xc_func_info_gga_x_ak13 = {
  XC_GGA_X_AK13,
  XC_EXCHANGE,
  "Armiento & Kuemmel 2013",
  XC_FAMILY_GGA,
  {&xc_ref_Armiento2013_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-24,
  2, ext_params, set_ext_params,
  gga_x_ak13_init, NULL,
  NULL, work_gga, NULL
};
