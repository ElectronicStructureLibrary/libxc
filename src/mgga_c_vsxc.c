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

#define XC_MGGA_C_VSXC          232 /* VSxc from Van Voorhis and Scuseria (correlation part) */

typedef struct{
  const double alpha_ss, alpha_ab;
  const double dss[6], dab[6];
} mgga_c_vsxc_params;

static const mgga_c_vsxc_params par_vsxc = {
  0.00515088, 0.00304966,
  { 3.270912e-01, -3.228915e-02, -2.942406e-02,  2.134222e-03, -5.451559e-03,  1.577575e-02},
  { 7.035010e-01,  7.694574e-03,  5.152765e-02,  3.394308e-05, -1.269420e-03,  1.296118e-03}
};

static void 
mgga_c_vsxc_init(xc_func_type *p)
{
  mgga_c_vsxc_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_vsxc_params));
  params = (mgga_c_vsxc_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_C_VSXC:
    memcpy(params, &par_vsxc, sizeof(mgga_c_vsxc_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_vsxc\n");
    exit(1);
  }  
}

#include "maple2c/mgga_c_vsxc.c"

#define func maple2c_func
#include "work_mgga_c.c"

const xc_func_info_type xc_func_info_mgga_c_vsxc = {
  XC_MGGA_C_VSXC,
  XC_CORRELATION,
  "VSXC (correlation part)",
  XC_FAMILY_MGGA,
  {&xc_ref_VanVoorhis1998_400, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1.0e-23,
  0, NULL, NULL,
  mgga_c_vsxc_init,
  NULL,
  NULL, NULL,
  work_mgga_c,
};
