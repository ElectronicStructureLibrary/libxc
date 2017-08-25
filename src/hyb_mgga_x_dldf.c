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

#define XC_HYB_MGGA_X_DLDF      36 /* Dispersionless Density Functional */

static void
mgga_x_dldf_init(xc_func_type *p)
{
  p->cam_alpha   = 0.6144129;
}

#include "maple2c/hyb_mgga_x_dldf.c"

#define func maple2c_func
#include "work_mgga_x.c"

const xc_func_info_type xc_func_info_hyb_mgga_x_dldf = {
  XC_HYB_MGGA_X_DLDF,
  XC_EXCHANGE,
  "Dispersionless Density Functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Pernal2009_263201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  MIN_DENS,
  0, NULL, NULL,
  mgga_x_dldf_init, NULL,
  NULL, NULL, work_mgga_x,
};
