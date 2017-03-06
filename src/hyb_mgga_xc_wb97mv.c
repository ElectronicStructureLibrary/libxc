/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

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

#define XC_HYB_MGGA_XC_WB97M_V   531 /* Mardirossian and Head-Gordon */

static void
hyb_mgga_xc_b97mv_init(XC(func_type) *p)
{
  p->cam_omega =  0.3;
  p->cam_alpha =  1.0;
  p->cam_beta  = -(1.0 - 0.15);
  p->nlc_b = 6.0;
  p->nlc_C = 0.01;
}

#include "maple2c/hyb_mgga_xc_b97mv.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_hyb_mgga_xc_wb97m_v) = {
  XC_HYB_MGGA_XC_WB97M_V,
  XC_EXCHANGE_CORRELATION,
  "wB97M-V exchange-correlation functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Mardirossian2016_214110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_b97mv_init,
  NULL, NULL, NULL,
  work_mgga_c,
};
