/*
 Copyright (C) 2015 Susi Lehtola

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

#define XC_GGA_XC_VV10         255 /* Vydrov and Van Voorhis */
#define XC_HYB_GGA_XC_LC_VV10  469 /* Vydrov and Van Voorhis */

static void
gga_xc_vv10_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_RPW86, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0, 1.0};
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  
  p->nlc_b = 5.9;
  p->nlc_C = 0.0093;
}

static void
hyb_gga_xc_lc_vv10_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0, 1.0};
  
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  
  p->cam_omega = 0.45;
  p->cam_alpha =  1.0;
  p->cam_beta  = -1.0;
  p->nlc_b = 6.3;
  p->nlc_C = 0.0089;
  XC(gga_x_wpbeh_set_params)(p->func_aux[0], p->cam_omega);  
}

const XC(func_info_type) XC(func_info_gga_xc_vv10) = {
  XC_GGA_XC_VV10,
  XC_EXCHANGE_CORRELATION,
  "Vydrov and Van Voorhis",
  XC_FAMILY_GGA,
  {&xc_ref_Vydrov2010_244103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  gga_xc_vv10_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_lc_vv10) = {
  XC_HYB_GGA_XC_LC_VV10,
  XC_EXCHANGE_CORRELATION,
  "Vydrov and Van Voorhis",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Vydrov2010_244103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_VV10 | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_lc_vv10_init,
  NULL, NULL, NULL, NULL
};
