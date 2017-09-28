/*
 Copyright (C) 2008 Georg Madsen

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

#define XC_MGGA_C_KCIS         562 /* Krieger, Chen, Iafrate, and Savin */
#define XC_HYB_MGGA_XC_B0KCIS  563 /* Hybrid based on KCIS */

#include "maple2c/mgga_c_kcis.c"

#define func maple2c_func
#include "work_mgga_c.c"

const xc_func_info_type xc_func_info_mgga_c_kcis = {
  XC_MGGA_C_KCIS,
  XC_CORRELATION,
  "Krieger, Chen, Iafrate, and Savin",
  XC_FAMILY_MGGA,
  {&xc_ref_Rey1998_581, &xc_ref_Krieger1999_463, &xc_ref_Krieger2001_48, &xc_ref_Toulouse2002_10465, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_DEVELOPMENT,
  1e-24,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, NULL, work_mgga_c
};

/*************************************************************/
void
xc_hyb_mgga_xc_b0kcis_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_MGGA_C_KCIS};
  static double funcs_coef[2] = {1.0 - 0.25, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

const xc_func_info_type xc_func_info_hyb_mgga_xc_b0kcis = {
  XC_HYB_MGGA_XC_B0KCIS,
  XC_EXCHANGE_CORRELATION,
  "Hybrid based on KCIS",
  XC_FAMILY_MGGA,
  {&xc_ref_Toulouse2002_10465, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_DEVELOPMENT,
  1e-32,
  0, NULL, NULL,
  xc_hyb_mgga_xc_b0kcis_init, NULL, 
  NULL, NULL, work_mgga_c
};
