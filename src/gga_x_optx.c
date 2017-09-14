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

#define XC_GGA_X_OPTX         110 /* Handy & Cohen OPTX 01                          */

typedef struct{
  double a, b, gamma;
} gga_x_optx_params;


static void 
gga_x_optx_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_optx_params));

  xc_gga_x_optx_set_params(p, 1.05151, 1.43169/X_FACTOR_C, 0.006);
}


void 
xc_gga_x_optx_set_params(xc_func_type *p, double a, double b, double gamma)
{
  gga_x_optx_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_optx_params *) (p->params);

  params->a     = a;
  params->b     = b;
  params->gamma = gamma;
}

#include "maple2c/gga_x_optx.c"

#define func maple2c_func
#include "work_gga_x.c"

const xc_func_info_type xc_func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  {&xc_ref_Handy2001_403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-22,
  0, NULL, NULL,
  gga_x_optx_init, NULL, 
  NULL, work_gga_x, NULL
};
