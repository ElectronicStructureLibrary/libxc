/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_AM05         120 /* Armiento & Mattsson 05 exchange                */

typedef struct{
  double alpha, c;
} gga_x_am05_params;

static void
gga_x_am05_init(xc_func_type *p)
{
  gga_x_am05_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_am05_params));
  params = (gga_x_am05_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/gga_exc/gga_x_am05.c"
#include "work_gga.c"

static const func_params_type ext_params[] = {
  {"_alpha", 2.804L, "alpha"},
  {"_c", 0.7168L, "c"}
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_am05_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_am05_params *) (p->params);

  params->alpha = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c     = get_ext_param(p->info->ext_params, ext_params, 1);
}

const xc_func_info_type xc_func_info_gga_x_am05 = {
  XC_GGA_X_AM05,
  XC_EXCHANGE,
  "Armiento & Mattsson 05",
  XC_FAMILY_GGA,
  {&xc_ref_Armiento2005_085108, &xc_ref_Mattsson2008_084714, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  2, ext_params, set_ext_params,
  gga_x_am05_init, NULL,
  NULL, work_gga, NULL
};
