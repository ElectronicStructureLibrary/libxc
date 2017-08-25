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

#define XC_HYB_MGGA_X_M05      438 /* M05 functional from Minnesota     */
#define XC_HYB_MGGA_X_M05_2X   439 /* M05-2X functional from Minnesota  */
#define XC_HYB_MGGA_X_M06_2X   450 /* M06-2X functional from Minnesota  */

typedef struct{
  double csi_HF;
  const double a[12];
} mgga_x_m05_params;

static const mgga_x_m05_params par_m05 = {
  1.0 - 0.28,
  {1.0, 0.08151, -0.43956, -3.22422, 2.01819, 8.79431, -0.00295,
   9.82029, -4.82351, -48.17574, 3.64802, 34.02248}
};

static const mgga_x_m05_params par_m05_2x = {
  1.0 - 0.56,
  {1.0, -0.56833, -1.30057, 5.50070, 9.06402, -32.21075, -23.73298,
   70.22996, 29.88614, -60.25778, -13.22205, 15.23694}
};

static const mgga_x_m05_params par_m06_2x = {
  1.0, /* the mixing is already included in the params->a */
  {4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00, -1.422133e+01,
   1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,  1.522808e+01,  1.522227e+01}
};

static void
mgga_x_m05_init(XC(func_type) *p)
{
  mgga_x_m05_params *params;

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m05_params));
  params = (mgga_x_m05_params *) (p->params);

  switch(p->info->number){
  case XC_HYB_MGGA_X_M05: 
    memcpy(params, &par_m05, sizeof(mgga_x_m05_params));
    p->cam_alpha   = 0.28;
    break;
  case XC_HYB_MGGA_X_M05_2X:
    memcpy(params, &par_m05_2x, sizeof(mgga_x_m05_params));
    p->cam_alpha   = 0.56;
    break;
  case XC_HYB_MGGA_X_M06_2X:
    memcpy(params, &par_m06_2x, sizeof(mgga_x_m05_params));
    p->cam_alpha   = 0.54;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_m05\n");
    exit(1);
  }

}

#include "maple2c/hyb_mgga_x_m05.c"

#define func maple2c_func
#include "work_mgga_x.c"


const XC(func_info_type) XC(func_info_hyb_mgga_x_m05) = {
  XC_HYB_MGGA_X_M05,
  XC_EXCHANGE,
  "Minnesota M05 functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_161103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_m05_init, NULL, 
  NULL, NULL, work_mgga_x,
};


const XC(func_info_type) XC(func_info_hyb_mgga_x_m05_2x) = {
  XC_HYB_MGGA_X_M05_2X,
  XC_EXCHANGE,
  "M05-2X functional from Minnesota",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2006_364, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_m05_init, NULL, 
  NULL, NULL, work_mgga_x,
};

const XC(func_info_type) XC(func_info_hyb_mgga_x_m06_2x) = {
  XC_HYB_MGGA_X_M06_2X,
  XC_EXCHANGE,
  "M06-2X functional from Minnesota",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  MIN_DENS,
  0, NULL, NULL,
  mgga_x_m05_init, NULL,
  NULL, NULL, work_mgga_x,
};
