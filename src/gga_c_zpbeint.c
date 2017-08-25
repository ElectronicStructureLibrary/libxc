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

#define XC_GGA_C_ZPBEINT       61 /* spin-dependent gradient correction to PBEint       */
#define XC_GGA_C_ZPBESOL       63 /* spin-dependent gradient correction to PBEsol       */

typedef struct{
  double beta, alpha;
} gga_c_zpbeint_params;

static void 
gga_c_zpbeint_init(XC(func_type) *p)
{
  gga_c_zpbeint_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_c_zpbeint_params));
  params = (gga_c_zpbeint_params *) (p->params);
 
  switch(p->info->number){
  case XC_GGA_C_ZPBEINT:
    params->beta  = 0.052;
    params->alpha = 2.4;
    break;
  case XC_GGA_C_ZPBESOL:
    params->beta  = 0.046;
    params->alpha = 4.8;
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_zpbeint\n");
    exit(1);
  }
}

#include "maple2c/gga_c_zpbeint.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_zpbeint) = {
  XC_GGA_C_ZPBEINT,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEint",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12,
  0, NULL, NULL,
  gga_c_zpbeint_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_zpbesol) = {
  XC_GGA_C_ZPBESOL,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEsol",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12,
  0, NULL, NULL,
  gga_c_zpbeint_init, NULL, 
  NULL, work_gga_c, NULL
};
