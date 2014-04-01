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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_HYB_GGA_XC_PBEH    406 /* aka PBE0 or PBE1PBE */
#define XC_HYB_GGA_XC_PBE0_13 456 /* PBE0-1/3            */

static void
hyb_gga_xc_pbeh_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

void 
XC(hyb_gga_xc_pbeh_set_params)(XC(func_type) *p, FLOAT alpha)
{
  assert(alpha>=0 && alpha<=1.0);

  p->cam_alpha   = alpha;
  p->mix_coef[0] = 1.0 - alpha;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbeh) = {
  XC_HYB_GGA_XC_PBEH,
  XC_EXCHANGE_CORRELATION,
  "PBEH (PBE0)",
  XC_FAMILY_HYB_GGA,
  "C Adamo and V Barone, J. Chem. Phys. 110, 6158 (1999)\n"
  "M. Ernzerhof, G. E. Scuseria, J. Chem. Phys. 110, 5029 (1999)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_pbeh_init,
  NULL, NULL, NULL, NULL /* this is taken care by the generic routine */
};


static void
hyb_gga_xc_pbe0_13_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0 - 1.0/3.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 1.0/3.0;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbe0_13) = {
  XC_HYB_GGA_XC_PBE0_13,
  XC_EXCHANGE_CORRELATION,
  "PBE0-1/3",
  XC_FAMILY_HYB_GGA,
  "P Cortona, J. Chem. Phys. 136, 086101 (2012)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_pbe0_13_init,
  NULL, NULL, NULL, NULL /* this is taken care by the generic routine */
};


