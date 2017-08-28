/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_X_TM          540 /* Tao and Mo 2016 */

#include "maple2c/mgga_x_tm.c"

#define func xc_mgga_x_tm_enhance
#include "work_mgga_x.c"

const xc_func_info_type xc_func_info_mgga_x_tm = {
  XC_MGGA_X_TM,
  XC_EXCHANGE,
  "Tao and Mo 2016",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2016_073001, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  5.0e-13,
  0, NULL, NULL,
  NULL, NULL,
  NULL, NULL, work_mgga_x,
};
