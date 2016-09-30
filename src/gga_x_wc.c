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

#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_WC         118 /* Wu & Cohen */

static FLOAT mu, c, kappa;

void gga_x_wc_init(XC(func_type) *p_)
{
  mu    = 0.2195149727645171;
  c     = (146.0/2025.0)*(4.0/9.0) - (73.0/405.0)*(2.0/3.0) + (mu - 10.0/81.0);
  kappa = 0.8040;
}

#include "hand_written/gga_x_wc.c"
#include "math2c/gga_x_wc.c"

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_wc) = {
  XC_GGA_X_WC,
  XC_EXCHANGE,
  "Wu & Cohen",
  XC_FAMILY_GGA,
  {&xc_ref_Wu2006_235116, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_wc_init, NULL, 
  NULL, work_gga_x, NULL
};

