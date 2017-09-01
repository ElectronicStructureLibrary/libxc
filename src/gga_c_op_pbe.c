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

#define XC_GGA_C_OP_PBE      86 /* one-parameter progressive functional (PBE version)     */

#include "maple2c/gga_c_op_pbe.c"

#define func maple2c_func
#include "work_gga_c.c"

const xc_func_info_type xc_func_info_gga_c_op_pbe = {
  XC_GGA_C_OP_PBE,
  XC_CORRELATION,
  "one-parameter progressive functional (PBE version)",
  XC_FAMILY_GGA,
  {&xc_ref_Tsuneda1999_10664, &xc_ref_Tsuneda1999_5656, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_c, NULL
};
