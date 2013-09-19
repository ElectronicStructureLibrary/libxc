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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define  XC_HYB_GGA_XC_CAMY_BLYP 455  /* BLYP with yukawa screening */

void
XC(hyb_gga_xc_camy_blyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_B88, XC_GGA_X_SFAT, XC_GGA_C_LYP};
  static FLOAT funcs_coef[3];

    p->cam_alpha = 0.20;
    p->cam_beta  = 0.80;
    p->cam_omega = 0.44;	/* we use omega for gamma here, 'cause
				   both denote dampening parameters for
				   range related interactions */

  funcs_coef[0] = 1.0 - p->cam_alpha - p->cam_beta;
  funcs_coef[1] = p->cam_beta; /* 1.0 - p->cam_alpha - p->cam_beta? */
  funcs_coef[2] = 1.0;

  XC(mix_init)(p, 3, funcs_id, funcs_coef);

  XC(gga_x_sfat_set_params)(p->func_aux[1], XC_GGA_X_B88, p->cam_omega);
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_camy_blyp) = {
  XC_HYB_GGA_XC_CAMY_BLYP,
  XC_EXCHANGE_CORRELATION,
  "CAMY version of BLYP",
  XC_FAMILY_HYB_GGA,
  "Y Akinaga, S Ten-no, Chem. Phys. Lett. 462, 348-351 (2008)",
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  XC(hyb_gga_xc_camy_blyp_init),
  NULL, NULL, NULL, NULL
};

