/*
 Copyright (C) 2006-2007 M.A.L. Marques and
                    2015 Susi Lehtola

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

#define XC_HYB_GGA_XC_EDF2        476 /* Empirical functional from Lin, George and Gill */

static void
hyb_gga_xc_edf2_init(XC(func_type) *p)
{
  static int   funcs_id  [6] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_LYP, XC_GGA_C_LYP};
  static FLOAT funcs_coef[6] = {0.2811, 0.6227, -0.0551, 0.3029, 0.5998, -0.0053};

  XC(mix_init)(p, 6, funcs_id, funcs_coef);  
  XC(gga_x_b88_set_params)(p->func_aux[2], 0.0035, 6.0);
  XC(gga_c_lyp_set_params)(p->func_aux[5], 0.055, 0.158, 0.25, 0.3505);
  p->cam_alpha = 0.1695;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_edf2) = {
  XC_HYB_GGA_XC_EDF2,
  XC_EXCHANGE_CORRELATION,
  "EDF2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Lin2004_365, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_edf2_init, 
  NULL, NULL, NULL, NULL
};
