/*
 Copyright (C) 2017 M.A.L. Marques

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

#define XC_GGA_C_ZVPBEINT       557 /* another spin-dependent correction to PBEint       */
#define XC_GGA_C_ZVPBESOL       558 /* another spin-dependent correction to PBEsol       */

typedef struct{
  FLOAT beta, alpha, omega;
} gga_c_zvpbeint_params;

static void 
gga_c_zvpbeint_init(XC(func_type) *p)
{
  gga_c_zvpbeint_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_c_zvpbeint_params));
  params = (gga_c_zvpbeint_params *) (p->params);
 
  switch(p->info->number){
  case XC_GGA_C_ZVPBEINT:
    params->beta  = 0.052;
    params->alpha = 1.0;
    params->omega = 4.5;
    break;
  case XC_GGA_C_ZVPBESOL:
    params->beta  = 0.046;
    params->alpha = 1.8;
    params->omega = 4.5;
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_zvpbeint\n");
    exit(1);
  }
}

#include "maple2c/gga_c_zvpbeint.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_zvpbeint) = {
  XC_GGA_C_ZVPBEINT,
  XC_CORRELATION,
  "another spin-dependent correction to PBEint",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2012_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12,
  0, NULL, NULL,
  gga_c_zvpbeint_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_zvpbesol) = {
  XC_GGA_C_ZVPBESOL,
  XC_CORRELATION,
  "another spin-dependent correction to PBEsol",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2012_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12,
  0, NULL, NULL,
  gga_c_zvpbeint_init, NULL, 
  NULL, work_gga_c, NULL
};
