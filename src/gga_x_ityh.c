/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_ITYH 529 /* short-range recipe B88 functionals - erf */

static const func_params_type ext_params[] = {
  {"_omega", 0.2, "Screening parameter"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  p->cam_omega = get_ext_param(p->info->ext_params, ext_params, 0);
}

#include "maple2c/gga_exc/gga_x_ityh.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_ityh = {
  XC_GGA_X_ITYH,
  XC_EXCHANGE,
  "Short-range recipe for B88 functional - erf",
  XC_FAMILY_GGA,
  {&xc_ref_Iikura2001_3540, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-8,
  1, ext_params, set_ext_params,
  NULL, NULL, 
  NULL, work_gga, NULL
};
