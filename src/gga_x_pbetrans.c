/*
 Copyright (C) 2017 Susi Lehtola

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

#define XC_GGA_X_PBETRANS  291 /* Gradient-based interpolation between PBE and revPBE */

#include "maple2c/gga_x_pbetrans.c"

#define func XC(gga_x_pbetrans_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pbetrans) = {
  XC_GGA_X_PBETRANS,
  XC_EXCHANGE,
  "Gradient-regulated connection-based correction for the PBE exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Bremond2016_1059, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL
};
