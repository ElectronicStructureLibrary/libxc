/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_OPTC       200 /* Optimized correlation functional of Cohen and Handy */

typedef struct{
  double c1, c2;
} gga_c_optc_params;

static void
gga_c_optc_init(xc_func_type *p)
{
  gga_c_optc_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_c_optc_params));
  params = (gga_c_optc_params *) (p->params);

  /* defaults set by set_ext_params */
}

#include "maple2c/gga_exc/gga_c_optc.c"
#include "work_gga.c"

static const func_params_type ext_params[] = {
  {"_c1", 1.1015L, "c1"},
  {"_c2", 0.6625L, "c2"}
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_c_optc_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_c_optc_params *) (p->params);

  params->c1 = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c2 = get_ext_param(p->info->ext_params, ext_params, 1);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_optc = {
  XC_GGA_C_OPTC,
  XC_CORRELATION,
  "Optimized correlation functional of Cohen and Handy",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2001_607, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-12,
  2, ext_params, set_ext_params,
  gga_c_optc_init, NULL,
  NULL, work_gga, NULL
};
