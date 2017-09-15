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

#define XC_MGGA_XC_ZLP          42 /* Zhao, Levy & Parr, Eq. (21) */

#include "maple2c/mgga_xc_zlp.c"

#define func maple2c_func
#include "work_mgga_c.c"

const xc_func_info_type xc_func_info_mgga_xc_zlp = {
  XC_MGGA_XC_ZLP,
  XC_EXCHANGE_CORRELATION,
  "Zhao, Levy & Parr, Eq. (21)",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao1993_918, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-23,
  0, NULL, NULL,
  NULL,
  NULL, NULL, NULL,
  work_mgga_c,
};

