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

#define XC_GGA_X_KT1          145 /* Exchange part of Keal and Tozer version 1 */
#define XC_GGA_XC_KT1         167 /* Keal and Tozer version 1                  */
#define XC_GGA_XC_KT2         146 /* Keal and Tozer version 2                  */

typedef struct{
  double gamma, delta;
} gga_x_kt_params;

static void 
gga_x_kt_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_kt_params));

  xc_gga_x_kt_set_params(p, -0.006, 0.1);
}


void 
xc_gga_x_kt_set_params(xc_func_type *p, double gamma, double delta)
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

const xc_func_info_type xc_func_info_gga_x_kt1 = {
  XC_GGA_X_KT1,
  XC_EXCHANGE,
  "Exchange part of Keal and Tozer version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_kt_init, NULL, 
  NULL, work_gga_c, NULL
};


static void
gga_xc_kt1_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);  
}

const xc_func_info_type xc_func_info_gga_xc_kt1 = {
  XC_GGA_XC_KT1,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 1",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_xc_kt1_init, NULL, 
  NULL, NULL, NULL
};

static void
gga_xc_kt2_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_KT1, XC_LDA_C_VWN};
  static double funcs_coef[3] = {1.07173 - 1.0, 1.0, 0.576727};

  xc_mix_init(p, 3, funcs_id, funcs_coef);  
}

const xc_func_info_type xc_func_info_gga_xc_kt2 = {
  XC_GGA_XC_KT2,
  XC_EXCHANGE_CORRELATION,
  "Keal and Tozer, version 2",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2003_3015, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_xc_kt2_init, NULL, 
  NULL, NULL, NULL
};
