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

#define XC_GGA_C_LYP    131  /* Lee, Yang & Parr */
#define XC_GGA_C_TM_LYP 559  /* Takkar and McCarthy reparametrization */

typedef struct{
  FLOAT A, B, c, d;
} gga_c_lyp_params;

void XC(gga_c_lyp_init)(XC(func_type) *p)
{
  assert(p->params == NULL);

  p->params = malloc(sizeof(gga_c_lyp_params));

  /* values of constants in standard LYP functional */
  switch(p->info->number){
  case XC_GGA_C_LYP:
    XC(gga_c_lyp_set_params)(p, 0.04918, 0.132, 0.2533, 0.349);
    break;
  case XC_GGA_C_TM_LYP:
    XC(gga_c_lyp_set_params)(p, 0.0393, 0.21, 0.41, 0.15);
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_pbe\n");
    exit(1);
  }
}


void XC(gga_c_lyp_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d)
{
  gga_c_lyp_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_c_lyp_params *) (p->params);

  params->A = A;
  params->B = B;
  params->c = c;
  params->d = d;
}


#include "maple2c/gga_c_lyp.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_lyp) = {
  XC_GGA_C_LYP,
  XC_CORRELATION,
  "Lee, Yang & Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Lee1988_785, &xc_ref_Miehlich1989_200, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  XC(gga_c_lyp_init), NULL,
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_tm_lyp) = {
  XC_GGA_C_TM_LYP,
  XC_CORRELATION,
  "Takkar and McCarthy reparametrization",
  XC_FAMILY_GGA,
  {&xc_ref_Thakkar2009_134109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  XC(gga_c_lyp_init), NULL,
  NULL, work_gga_c, NULL
};
