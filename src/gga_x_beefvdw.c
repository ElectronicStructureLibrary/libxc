/*
 Copyright (C) 2014 Jess Wellendorff, M.A.L. Marques

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

#define XC_GGA_X_BEEFVDW          285 /* BEEF-vdW exchange */
#define XC_GGA_XC_BEEFVDW         286 /* BEEF-vdW exchange-correlation */

#include "maple2c/gga_x_beefvdw.c"

#define func maple2c_func
#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_beefvdw) = {
  XC_GGA_X_BEEFVDW,
  XC_EXCHANGE,
  "BEEF-vdW exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL,
};


void
gga_xc_beefvdw_init(XC(func_type) *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_BEEFVDW, XC_LDA_C_PW_MOD, XC_GGA_C_PBE};
  static FLOAT funcs_coef[3] = {1.0, 0.6001664769, 1.0 - 0.6001664769};

  XC(mix_init)(p, 3, funcs_id, funcs_coef);
}

const XC(func_info_type) XC(func_info_gga_xc_beefvdw) = {
  XC_GGA_XC_BEEFVDW,
  XC_EXCHANGE_CORRELATION,
  "BEEF-vdW exchange-correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Wellendorff2012_235149, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  gga_xc_beefvdw_init, NULL, 
  NULL, NULL, NULL,
};
