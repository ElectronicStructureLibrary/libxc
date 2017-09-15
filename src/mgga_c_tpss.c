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

#define XC_MGGA_C_TPSS          231 /* Tao, Perdew, Staroverov & Scuseria correlation */

typedef struct{
  double beta, d;
  double C0_c[4];
} mgga_c_tpss_params;


static void 
mgga_c_tpss_init(xc_func_type *p)
{

  assert(p != NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_tpss_params));

  switch(p->info->number){
  case XC_MGGA_C_TPSS:
    xc_mgga_c_tpss_set_params(p, 0.06672455060314922, 2.8, 0.53, 0.87, 0.50, 2.26);
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_tpss\n");
    exit(1);
  }
}

void
xc_mgga_c_tpss_set_params
     (xc_func_type *p, double beta, double d, double C0_0, double C0_1, double C0_2, double C0_3)
{
  mgga_c_tpss_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_c_tpss_params *) (p->params);

  params->beta    = beta;
  params->d       = d;
  params->C0_c[0] = C0_0;
  params->C0_c[1] = C0_1;
  params->C0_c[2] = C0_2;
  params->C0_c[3] = C0_3;
}

#include "maple2c/mgga_c_tpss.c"

#define func maple2c_func
#include "work_mgga_c.c"


const xc_func_info_type xc_func_info_mgga_c_tpss = {
  XC_MGGA_C_TPSS,
  XC_CORRELATION,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Tao2003_146401, &xc_ref_Perdew2004_6898, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-23, /* densities smaller than 1e-26 give NaNs */
  0, NULL, NULL,
  mgga_c_tpss_init,
  NULL, NULL, NULL,
  work_mgga_c,
};
