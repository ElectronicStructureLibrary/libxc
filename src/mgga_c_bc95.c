/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_C_BC95          240 /* Becke correlation 95 */

typedef struct{
  double css, copp;
} mgga_c_bc95_params;


static void 
mgga_c_bc95_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_bc95_params));

  xc_mgga_c_bc95_set_params(p, 0.038, 0.0031);
}


void 
xc_mgga_c_bc95_set_params(xc_func_type *p, double css, double copp)
{
  mgga_c_bc95_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_c_bc95_params *) (p->params);

  params->css  = css;
  params->copp = copp;
}


#include "maple2c/mgga_c_bc95.c"

#define func maple2c_func
#include "work_mgga_c.c"

const xc_func_info_type xc_func_info_mgga_c_bc95 = {
  XC_MGGA_C_BC95,
  XC_CORRELATION,
  "Becke correlation 95",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke1996_1040, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  mgga_c_bc95_init,
  NULL, NULL, NULL,
  work_mgga_c,
};

