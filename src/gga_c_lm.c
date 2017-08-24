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

/************************************************************************
  This functional is provided for historical reasons.
  It was one of the first GGAs that ever appeared.
************************************************************************/

#define XC_GGA_C_LM          137 /* Langreth and Mehl correlation          */

#include "maple2c/gga_c_lm.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_lm) = {
  XC_GGA_C_LM,
  XC_CORRELATION,
  "Langreth & Mehl",
  XC_FAMILY_GGA,
  {&xc_ref_Langreth1981_446, &xc_ref_Hu1985_391, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-10, 1e-32, 0.0,
  0, NULL, NULL,
  NULL, NULL,
  NULL, work_gga_c, NULL
};

