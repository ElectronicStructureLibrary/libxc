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

#define XC_MGGA_X_MS0          221 /* MS exchange of Sun, Xiao, and Ruzsinszky */
#define XC_MGGA_X_MS1          222 /* MS1 exchange of Sun, et al */
#define XC_MGGA_X_MS2          223 /* MS2 exchange of Sun, et al */
#define XC_HYB_MGGA_X_MS2H     224 /* MS2 hybrid exchange of Sun, et al */

typedef struct{
  double kappa, c, b;
} mgga_x_ms_params;

static void 
mgga_x_ms_init(XC(func_type) *p)
{
  mgga_x_ms_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_ms_params));
  params = (mgga_x_ms_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_MS0:
    params->kappa = 0.29;
    params->c     = 0.28771;
    params->b     = 1.0;
    break;
  case XC_MGGA_X_MS1:
    params->kappa = 0.404;
    params->c     = 0.18150;
    params->b     = 1.0;
    break;
  case XC_MGGA_X_MS2:
    params->kappa = 0.504;
    params->c     = 0.14601;
    params->b     = 4.0;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_ms\n");
    exit(1);
  }
}

#include "maple2c/mgga_x_ms.c"

#define func maple2c_func
#include "work_mgga_x.c"


const XC(func_info_type) XC(func_info_mgga_x_ms0) = {
  XC_MGGA_X_MS0,
  XC_EXCHANGE,
  "MS exchange of Sun, Xiao, and Ruzsinszky",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2012_051101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_ms_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_ms1) = {
  XC_MGGA_X_MS1,
  XC_EXCHANGE,
  "MS1 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_ms_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_ms2) = {
  XC_MGGA_X_MS2,
  XC_EXCHANGE,
  "MS2 exchange of Sun, et al",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_ms_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

static void
hyb_mgga_x_ms2h_init(XC(func_type) *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_MS2};
  static double funcs_coef[1] = {0.91};

  XC(mix_init)(p, 1, funcs_id, funcs_coef);
  p->cam_alpha = 0.09;
}


const XC(func_info_type) XC(func_info_hyb_mgga_x_ms2h) = {
  XC_HYB_MGGA_X_MS2H,
  XC_EXCHANGE,
  "MS2 hybrid exchange of Sun, et al",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Sun2013_044113, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_x_ms2h_init, NULL, 
  NULL, NULL, NULL
};
