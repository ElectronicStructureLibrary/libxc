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

#define XC_MGGA_C_TPSSLOC       247 /* Semilocal dynamical correlation */

#include "maple2c/mgga_c_tpssloc.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_tpssloc) = {
  XC_MGGA_C_TPSSLOC,
  XC_CORRELATION,
  "Semilocal dynamical correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Constantin2012_035130, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-26, 1e-32, /* densities smaller than 1e-26 give NaNs */
  0, NULL, NULL,
  NULL, NULL,
  NULL, NULL, work_mgga_c
};

