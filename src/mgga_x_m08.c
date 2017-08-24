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

#define XC_MGGA_X_M08_HX       219 /* Worker for M08-HX functional     */
#define XC_MGGA_X_M08_SO       220 /* Worker for M08-SO functional     */
#define XC_HYB_MGGA_XC_M08_HX  460 /* M08-HX functional from Minnesota */
#define XC_HYB_MGGA_XC_M08_SO  461 /* M08-SO functional from Minnesota */

typedef struct{
  const FLOAT a[12], b[12];
} mgga_x_m08_params;

static const mgga_x_m08_params par_m08_hx = {
  {
     1.3340172e+00, -9.4751087e+00, -1.2541893e+01,  9.1369974e+00,  3.4717204e+01,  5.8831807e+01,
     7.1369574e+01,  2.3312961e+01,  4.8314679e+00, -6.5044167e+00, -1.4058265e+01,  1.2880570e+01
  }, {
    -8.5631823e-01,  9.2810354e+00,  1.2260749e+01, -5.5189665e+00, -3.5534989e+01, -8.2049996e+01,
    -6.8586558e+01,  3.6085694e+01, -9.3740983e+00, -5.9731688e+01,  1.6587868e+01,  1.3993203e+01
  }
};

static const mgga_x_m08_params par_m08_so = {
  {
    -3.4888428e-01, -5.8157416e+00,  3.7550810e+01,  6.3727406e+01, -5.3742313e+01, -9.8595529e+01,
     1.6282216e+01,  1.7513468e+01, -6.7627553e+00,  1.1106658e+01,  1.5663545e+00,  8.7603470e+00
  }, {
     7.8098428e-01,  5.4538178e+00, -3.7853348e+01, -6.2295080e+01,  4.6713254e+01,  8.7321376e+01,
     1.6053446e+01,  2.0126920e+01, -4.0343695e+01, -5.8577565e+01,  2.0890272e+01,  1.0946903e+01
  }
};


static void
mgga_x_m08_init(XC(func_type) *p)
{
  mgga_x_m08_params *params;

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m08_params));
  params = (mgga_x_m08_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_M08_HX: 
    memcpy(params, &par_m08_hx, sizeof(mgga_x_m08_params));
    break;
  case XC_MGGA_X_M08_SO:
    memcpy(params, &par_m08_so, sizeof(mgga_x_m08_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_m08\n");
    exit(1);
  }
}

#include "maple2c/mgga_x_m08.c"

#define func maple2c_func
#include "work_mgga_x.c"


const XC(func_info_type) XC(func_info_mgga_x_m08_hx) = {
  XC_MGGA_X_M08_HX,
  XC_EXCHANGE,
  "Worker for hyb_mgga_x_m08_hx",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_1849, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_m08_so) = {
  XC_MGGA_X_M08_SO,
  XC_EXCHANGE,
  "Worker for hyb_mgga_x_m08_so",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_1849, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

static void
hyb_mgga_xc_m08_hx_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_M08_HX, XC_MGGA_C_M08_HX};
  static FLOAT funcs_coef[2] = {1.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.5223;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_m08_hx) = {
  XC_HYB_MGGA_XC_M08_HX,
  XC_EXCHANGE_CORRELATION,
  "Minnesota M08-HX hybrid functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2008_1849, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_m08_hx_init,
  NULL, NULL, NULL, NULL
};

static void
hyb_mgga_xc_m08_so_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_M08_SO, XC_MGGA_C_M08_SO};
  static FLOAT funcs_coef[2] = {1.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.5679;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_m08_so) = {
  XC_HYB_MGGA_XC_M08_SO,
  XC_EXCHANGE_CORRELATION,
  "Minnesota M08-SO hybrid functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2008_1849, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_m08_so_init,
  NULL, NULL, NULL, NULL
};
