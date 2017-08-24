/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include "util.h"

#define XC_MGGA_X_M11          225 /* Worker for M11 functional        */
#define XC_HYB_MGGA_XC_M11     462 /* M11    functional from Minnesota */

typedef struct{
  const FLOAT a[12], b[12];
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


static void
mgga_x_m11_init(XC(func_type) *p)
{
  mgga_x_m11_params *params;

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m11_params));
  params = (mgga_x_m11_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_M11:
    memcpy(params, &par_m11, sizeof(mgga_x_m11_params));
    p->cam_omega = 0.25;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_m08\n");
    exit(1);
  }
}

#include "maple2c/mgga_x_m11.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_x_m11) = {
  XC_MGGA_X_M11,
  XC_EXCHANGE,
  "Worker for hyb_mgga_xc_m11",
  XC_FAMILY_MGGA,
  {&xc_ref_Peverati2011_2810, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  mgga_x_m11_init, NULL, 
  NULL, NULL, work_mgga_c,
};


static void
hyb_mgga_xc_m11_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_M11, XC_MGGA_C_M11};
  static FLOAT funcs_coef[2] = {1.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 1.0;
  p->cam_beta  = -(1.0 - 0.428);
  p->cam_omega = 0.25;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_m11) = {
  XC_HYB_MGGA_XC_M11,
  XC_EXCHANGE_CORRELATION,
  "Minnesota M11 hybrid functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Peverati2011_2810, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_m11_init,
  NULL, NULL, NULL, NULL
};


