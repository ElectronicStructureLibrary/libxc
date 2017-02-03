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

#define XC_GGA_X_KT1          145 /* Keal and Tozer version 1             */
#define XC_GGA_XC_KT2         146 /* Keal and Tozer version 2             */

typedef struct{
  FLOAT gamma, delta;
} gga_x_kt_params;

static void 
gga_x_kt_init(XC(func_type) *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_kt_params));

  XC(gga_x_kt_set_params)(p, -0.006, 0.1);
}


void 
XC(gga_x_kt_set_params)(XC(func_type) *p, FLOAT gamma, FLOAT delta)
{
  gga_x_kt_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_kt_params *) (p->params);

  params->gamma = gamma;
  params->delta = delta;
}

#include "maple2c/gga_x_kt.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_x_kt1) = {
  XC_GGA_X_KT1,
  XC_EXCHANGE,
  "Keal and Tozer, version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_kt_init, NULL, 
  NULL, work_gga_c, NULL
};


static void
gga_xc_kt2_init(XC(func_type) *p)
{
  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN};
  static FLOAT funcs_coef[3] = {1.07173 - 1.0, 1.0, 0.576727};

  XC(mix_init)(p, 3, funcs_id, funcs_coef);  
}

const XC(func_info_type) XC(func_info_gga_xc_kt2) = {
  XC_GGA_XC_KT2,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 2",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_xc_kt2_init, NULL, 
  NULL, NULL, NULL
};


