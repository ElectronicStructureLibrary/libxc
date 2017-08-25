/*
 Copyright (C) 2013 Rolf Wuerdemann, M.A.L. Marques

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

#define XC_HYB_GGA_XC_CAMY_B3LYP        470 /* B3LYP with Yukawa screening */

void
xc_hyb_gga_xc_camy_b3lyp_init(xc_func_type *p)
{
  static double ac = 0.81;
  static int   funcs_id  [4] = {XC_GGA_X_B88, XC_GGA_X_SFAT, XC_LDA_C_VWN, XC_GGA_C_LYP};
  static double funcs_coef[4];

  /* Need temp variables since cam_ parameters are initialized in mix_init */
  static double omega, alpha, beta;

  /* N.B. The notation used in the original reference uses a different
     convention for alpha and beta.  In libxc, alpha is the weight for
     HF exchange, which in the original reference is alpha+beta.
  */
  omega = 0.34;
  alpha = 0.65;
  beta  = -0.46;

  funcs_coef[0] = 1.0 - alpha;
  funcs_coef[1] = -beta;
  funcs_coef[2] = 1.0 - ac;
  funcs_coef[3] = ac;
  
  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_gga_x_sfat_set_params(p->func_aux[1], XC_GGA_X_B88, omega);
  
  p->cam_omega = omega;
  p->cam_alpha = alpha;
  p->cam_beta  = beta;
}

const xc_func_info_type xc_func_info_hyb_gga_xc_camy_b3lyp = {
  XC_HYB_GGA_XC_CAMY_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAMY version of B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Seth2012_901, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAMY | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  xc_hyb_gga_xc_camy_b3lyp_init,
  NULL, NULL, NULL, NULL
};

