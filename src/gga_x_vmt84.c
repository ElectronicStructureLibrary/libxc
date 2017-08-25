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

#define XC_GGA_X_VMT84_PBE        69 /* VMT{8,4} with constraint satisfaction with mu = mu_PBE  */
#define XC_GGA_X_VMT84_GE         68 /* VMT{8,4} with constraint satisfaction with mu = mu_GE  */

typedef struct{
  double mu;
  double alpha;
} gga_x_vmt84_params;

static void 
gga_x_vmt84_init(XC(func_type) *p)
{
  gga_x_vmt84_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_vmt84_params));
  params = (gga_x_vmt84_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_VMT84_PBE:
    params->mu    = 0.2195149727645171;
    params->alpha = 0.000074;
    break;
  case XC_GGA_X_VMT84_GE:
    params->mu = 10.0/81.0;
    params->alpha = 0.000023;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_vmt84\n");
    exit(1);
  }
}

#include "maple2c/gga_x_vmt84.c"

#define func maple2c_func
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_vmt84_pbe) = {
  XC_GGA_X_VMT84_PBE,
  XC_EXCHANGE,
  "VMT{8,4} with constraint satisfaction with mu = mu_PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Vela2012_144115, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_vmt84_init, NULL, 
  NULL, work_gga_x, NULL
};

const XC(func_info_type) XC(func_info_gga_x_vmt84_ge) = {
  XC_GGA_X_VMT84_GE,
  XC_EXCHANGE,
  "VMT{8,4} with constraint satisfaction with mu = mu_GE",
  XC_FAMILY_GGA,
  {&xc_ref_Vela2012_144115, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_vmt84_init, NULL, 
  NULL, work_gga_x, NULL
};
