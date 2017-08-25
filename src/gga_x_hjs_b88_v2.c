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

#define XC_GGA_X_HJS_B88_V2   46 /* HJS screened exchange corrected B88 version */

typedef struct{
  double omega;
} gga_x_hjs_b88_v2_params;

static void
gga_x_hjs_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_hjs_b88_v2_params));

  /* we take 0.11 as the default for hjs_b88_v2 */
  xc_gga_x_hjs_b88_v2_set_params(p, 0.11);
}

void 
xc_gga_x_hjs_b88_v2_set_params(xc_func_type *p, double omega)
{
  gga_x_hjs_b88_v2_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_hjs_b88_v2_params *) (p->params);

  params->omega = omega;
}

#include "maple2c/gga_x_hjs_b88_v2.c"

#define func maple2c_func
#include "work_gga_c.c"

const xc_func_info_type xc_func_info_gga_x_hjs_b88_v2 = {
  XC_GGA_X_HJS_B88_V2,
  XC_EXCHANGE,
  "HJS screened exchange B88 corrected version",
  XC_FAMILY_GGA,
  {&xc_ref_Weintraub2009_754, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-6, /* densities smaller than 1e-6 yield NaNs */
  0, NULL, NULL,
  gga_x_hjs_init, NULL, 
  NULL, work_gga_c, NULL
};
