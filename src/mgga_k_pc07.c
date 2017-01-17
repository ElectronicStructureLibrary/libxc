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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_MGGA_K_PC07          543 /* Perdew and Constantin 2007 */

#include "maple2c/mgga_k_pc07.c"

#define func XC(mgga_k_pc07_enhance)
#define XC_KINETIC_FUNCTIONAL
#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_k_pc07) = {
  XC_MGGA_K_PC07,
  XC_EXCHANGE,
  "Perdew and Constantin 2007",
  XC_FAMILY_MGGA,
  {&xc_ref_Perdew2007_155109, NULL, NULL, NULL, NULL},
  XC_FLAGS_DEVELOPMENT | XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  0, NULL, NULL,
  NULL, NULL,
  NULL, NULL, work_mgga_k,
};
