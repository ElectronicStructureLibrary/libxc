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

#define  XC_HYB_GGA_XC_LCY_BLYP 468  /* BLYP with yukawa screening */

void
XC(hyb_gga_xc_lcy_blyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_SFAT, XC_GGA_C_LYP};
  static FLOAT funcs_coef[2];

  FLOAT gamma = 0.75;
  
  funcs_coef[0] = 1.0;
  funcs_coef[1] = 1.0;

  XC(mix_init)(p, 2, funcs_id, funcs_coef);

  XC(gga_x_sfat_set_params)(p->func_aux[0], XC_GGA_X_B88, gamma);
  p->cam_omega = gamma;
  p->cam_alpha = 1.0;
  p->cam_beta  = -1.0;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_lcy_blyp) = {
  XC_HYB_GGA_XC_LCY_BLYP,
  XC_EXCHANGE_CORRELATION,
  "LCY version of BLYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Akinaga2008_348, &xc_ref_Seth2013_2286, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_LCY | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  XC(hyb_gga_xc_lcy_blyp_init),
  NULL, NULL, NULL, NULL
};

