/*
 Copyright (C) 2013 Rolf Wuerdemann, M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_SFAT 530 /* short-range recipe for PBE functional */
/* see
   Savin, Flad, Int. J. Quant. Chem. Vol. 56, 327-332 (1995)
   Akinaga, Ten-no, Chem. Phys. Lett. 462 (2008) 348-351
*/

static const func_params_type ext_params[] = {
  {"_omega", 0.44, "Screening parameter"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  p->cam_omega = get_ext_param(p->info->ext_params, ext_params, 0);
}

#include "maple2c/gga_exc/gga_x_sfat.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_sfat = {
  XC_GGA_X_SFAT,
  XC_EXCHANGE,
  "Short-range recipe for B88 functional - Yukawa",
  XC_FAMILY_GGA,
  {&xc_ref_Savin1995_327, &xc_ref_Akinaga2008_348, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-19,
  1, ext_params, set_ext_params,
  NULL, NULL,
  NULL, work_gga, NULL
};

