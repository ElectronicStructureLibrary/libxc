/*
 Copyright (C) 2008 Georg Madsen

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

#define XC_GGA_X_AIRY  192 /* Constantin et al based on the Airy gas */
#define XC_GGA_X_LAG   193 /* Local Airy Gas */

static FLOAT 
  a1  =   0.041106, a2  =   2.626712, a3  =   0.092070, a4  =   0.657946, 
  a5  = 133.983631, a6  =   3.217063, a7  = 136.707378, a8  =   3.223476, 
  a9  =   2.675484, a10 =   3.473804;


#include "hand_written/gga_x_airy.c"

#include "math2c/gga_x_lag.c"
#undef math2c_func
#include "math2c/gga_x_airy.c"
#undef math2c_func

static inline void math2c_new_func
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  switch(p->info->number){
  case XC_GGA_X_AIRY:
    XC(math2c_gga_x_airy_enhance)(p, order, x, f, dfdx, d2fdx2, d3fdx3);
    break;
  case XC_GGA_X_LAG:
    XC(math2c_gga_x_lag_enhance)(p, order, x, f, dfdx, d2fdx2, d3fdx3);
    break;
  }
}
#define math2c_func math2c_new_func

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_airy) = {
  XC_GGA_X_AIRY,
  XC_EXCHANGE,
  "Constantin et al based on the Airy gas",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2009_035125, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL
};

const XC(func_info_type) XC(func_info_gga_x_lag) = {
  XC_GGA_X_LAG,
  XC_EXCHANGE,
  "Local Airy Gas",
  XC_FAMILY_GGA,
  {&xc_ref_Vitos2000_10046, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL
};
