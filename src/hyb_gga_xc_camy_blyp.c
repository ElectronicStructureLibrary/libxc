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

#define  XC_HYB_GGA_XC_CAMY_BLYP 455  /* BLYP with yukawa screening */

void
XC(hyb_gga_xc_camy_blyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_B88, XC_GGA_X_SFAT, XC_GGA_C_LYP};
  static double funcs_coef[3];

  /* N.B. The notation used in the original reference uses a different
     convention for alpha and beta.  In libxc, alpha is the weight for
     HF exchange, which in the original reference is alpha+beta.
  */
  static double alpha, beta, omega;
  
  alpha = 1.00;
  beta  =-0.80;
  omega = 0.44;	/* we use omega for gamma here, 'cause
		   both denote dampening parameters for
	       	   range related interactions */
  
  funcs_coef[0] = 1.0 - alpha;
  funcs_coef[1] =-beta;
  funcs_coef[2] = 1.0;

  XC(mix_init)(p, 3, funcs_id, funcs_coef);
  XC(gga_x_sfat_set_params)(p->func_aux[1], XC_GGA_X_B88, omega);

  p->cam_omega=omega;
  p->cam_alpha=alpha;
  p->cam_beta=beta;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_camy_blyp) = {
  XC_HYB_GGA_XC_CAMY_BLYP,
  XC_EXCHANGE_CORRELATION,
  "CAMY version of BLYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Akinaga2008_348, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAMY | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_camy_blyp_init),
  NULL, NULL, NULL, NULL
};

