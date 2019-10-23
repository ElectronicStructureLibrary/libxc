/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_X_M11       297 /* M11 hybrid exchange functional from Minnesota        */
#define XC_HYB_MGGA_X_REVM11    304 /* revM11 hybrid exchange functional from Minnesota     */

typedef struct{
  const double a[12], b[12];
} mgga_x_m11_params;

static const mgga_x_m11_params par_m11 = {
  {
    -0.18399900e+00, -1.39046703e+01,  1.18206837e+01,  3.10098465e+01, -5.19625696e+01,  1.55750312e+01,
    -6.94775730e+00, -1.58465014e+02, -1.48447565e+00,  5.51042124e+01, -1.34714184e+01,  0.00000000e+00
  }, {
     0.75599900e+00,  1.37137944e+01, -1.27998304e+01, -2.93428814e+01,  5.91075674e+01, -2.27604866e+01,
    -1.02769340e+01,  1.64752731e+02,  1.85349258e+01, -5.56825639e+01,  7.47980859e+00,  0.00000000e+00
  }
};

static const mgga_x_m11_params par_revm11 = {
  {
   -0.3288860885e+00, -8.3888150476e+00,  0.7123891057e+00,  3.6196212952e+00,  4.3941708207e+00,  5.0453345584e+00,
    7.8667061191e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00
  }, {
    1.1038860885e+00,  8.0476369587e+00, -0.7353624773e+00, -2.4735275550e+00, -4.7319060355e+00, -5.8502502096e+00,
   -7.5059975327e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00,  0.0000000000e+00
  }
};

static void
mgga_x_m11_init(xc_func_type *p)
{
  mgga_x_m11_params *params;

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m11_params));
  params = (mgga_x_m11_params *) (p->params);

  switch(p->info->number){
  case XC_HYB_MGGA_X_M11:
    memcpy(params, &par_m11, sizeof(mgga_x_m11_params));
    p->cam_alpha = 1.0;
    p->cam_beta  = -(1.0 - 0.428);
    p->cam_omega = 0.25;
    break;
  case XC_HYB_MGGA_X_REVM11:
    memcpy(params, &par_revm11, sizeof(mgga_x_m11_params));
    p->cam_alpha = 1.0;
    p->cam_beta  = -(1.0 - 0.225);
    p->cam_omega = 0.40;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_m11\n");
    exit(1);
  }
}

#include "maple2c/mgga_exc/mgga_x_m11.c"
#include "work_mgga.c"

const xc_func_info_type xc_func_info_hyb_mgga_x_m11 = {
  XC_HYB_MGGA_X_M11,
  XC_EXCHANGE,
  "Minnesota M11 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Peverati2011_2810, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-11,
  0, NULL, NULL,
  mgga_x_m11_init, NULL,
  NULL, NULL, work_mgga,
};

const xc_func_info_type xc_func_info_hyb_mgga_x_revm11 = {
  XC_HYB_MGGA_X_REVM11,
  XC_EXCHANGE,
  "Revised Minnesota M11 hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Verma2019_2966, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-11,
  0, NULL, NULL,
  mgga_x_m11_init, NULL,
  NULL, NULL, work_mgga,
};
