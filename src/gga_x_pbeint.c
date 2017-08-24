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

#define XC_GGA_X_PBEINT        60 /* PBE for hybrid interfaces                      */
#define XC_GGA_K_APBEINT       54 /* interpolated version of APBE                   */
#define XC_GGA_K_REVAPBEINT    53 /* interpolated version of REVAPBE                */


typedef struct{
  FLOAT kappa, alpha, muPBE, muGE;
} gga_x_pbeint_params;


static void 
gga_x_pbe_init(XC(func_type) *p)
{
  gga_x_pbeint_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_pbeint_params));
  params = (gga_x_pbeint_params *) (p->params);
 
  switch(p->info->number){
  case XC_GGA_X_PBEINT:
    XC(gga_x_pbeint_set_params)(p, 0.8040, 0.197, 0.2195149727645171, MU_GE);
    break;
  case XC_GGA_K_APBEINT:
    XC(gga_x_pbeint_set_params)(p, 0.8040, 5.0/3.0, 0.23899, 5.0/27.0);
    break;
  case XC_GGA_K_REVAPBEINT:
    XC(gga_x_pbeint_set_params)(p, 1.245, 5.0/3.0, 0.23899, 5.0/27.0);
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_pbeint\n");
    exit(1);
  }
}


void 
XC(gga_x_pbeint_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT alpha, FLOAT muPBE, FLOAT muGE)
{
  gga_x_pbeint_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_pbeint_params *) (p->params);

  params->kappa = kappa;
  params->alpha = alpha;
  params->muPBE = muPBE;
  params->muGE  = muGE;
}

#include "maple2c/gga_x_pbeint.c"

#define func XC(gga_x_pbeint_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pbeint) = {
  XC_GGA_X_PBEINT,
  XC_EXCHANGE,
  "PBE for hybrid interfaces",
  XC_FAMILY_GGA,
  {&xc_ref_Fabiano2010_113104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga_x, NULL
};

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_apbeint) = {
  XC_GGA_K_APBEINT,
  XC_KINETIC,
  "interpolated version of APBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga_k, NULL
};

const XC(func_info_type) XC(func_info_gga_k_revapbeint) = {
  XC_GGA_K_REVAPBEINT,
  XC_KINETIC,
  "interpolated version of revAPBE",
  XC_FAMILY_GGA,
  {&xc_ref_Laricchia2011_2439, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  gga_x_pbe_init, NULL, 
  NULL, work_gga_k, NULL
};
