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

#include "util.h"

#define XC_GGA_X_2D_B86_MGC      124 /* Becke 86 MGC for 2D systems */

#include "maple2c/gga_x_2d_b86_mgc.c"

#define func maple2c_func
#define XC_DIMENSIONS 2
#include "work_gga_x.c"

const xc_func_info_type xc_func_info_gga_x_2d_b86_mgc = {
  XC_GGA_X_2D_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction for 2D",
  XC_FAMILY_GGA,
  {&xc_ref_Pittalis2009_012503, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-23,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL
};
