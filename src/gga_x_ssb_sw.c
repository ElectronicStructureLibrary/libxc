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

#define XC_GGA_X_SSB_SW       90  /* Swart, Sola and Bickelhaupt correction to PBE  */
#define XC_GGA_X_SSB          91  /* Swart, Sola and Bickelhaupt  */
#define XC_GGA_X_SSB_D        92  /* Swart, Sola and Bickelhaupt dispersion  */

typedef struct{
  FLOAT A, B, C, D, E;
} gga_x_ssb_sw_params;


static void 
gga_x_ssb_sw_init(XC(func_type) *p)
{

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ssb_sw_params));

  XC(gga_x_ssb_sw_set_params)(p, 1.0515, 0.191458, 0.254443, 0.180708, 4.036674);
}


void 
XC(gga_x_ssb_sw_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT C, FLOAT D, FLOAT E)
{
  gga_x_ssb_sw_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ssb_sw_params *) (p->params);

  params->A = A;
  params->B = B;
  params->C = C;
  params->D = D;
  params->E = E;
}

#include "maple2c/gga_x_ssb_sw.c"

#define func maple2c_func
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_ssb_sw) = {
  XC_GGA_X_SSB_SW,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt correction to PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_69, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_ssb_sw_init, NULL, 
  NULL, work_gga_x, NULL
};

static void
gga_x_ssb_init(XC(func_type) *p)
{
  static const FLOAT u = -1.205643, F = 0.995010, B = 0.137574;

  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_SSB_SW, XC_GGA_X_KT1};
  static FLOAT funcs_coef[3] = {-1.0, 1.0, 1.0};

  XC(mix_init)(p, 3, funcs_id, funcs_coef);  

  XC(gga_x_ssb_sw_set_params)(p->func_aux[1], 1.071769, 0.137574, 0.187883, 0.137574*(1.0 + 1.205643), 6.635315);
  XC(gga_x_kt_set_params)(p->func_aux[2], u*F*X_FACTOR_C*B*(X2S*X2S), 0.1);
}


const XC(func_info_type) XC(func_info_gga_x_ssb) = {
  XC_GGA_X_SSB,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_ssb_init, NULL, 
  NULL, NULL, NULL
};


static void
gga_x_ssb_d_init(XC(func_type) *p)
{
  static const FLOAT u = -0.749940, F = 0.949488, B = 0.197465;

  static int   funcs_id  [3] = {XC_LDA_X, XC_GGA_X_SSB_SW, XC_GGA_X_KT1};
  static FLOAT funcs_coef[3] = {-1.0, 1.0, 1.0};

  XC(mix_init)(p, 3, funcs_id, funcs_coef);  

  XC(gga_x_ssb_sw_set_params)(p->func_aux[1], 1.079966, 0.197465, 0.272729, 0.197465*(1.0 + 0.749940), 5.873645);
  XC(gga_x_kt_set_params)(p->func_aux[2], u*F*X_FACTOR_C*B*(X2S*X2S), 0.1);
}

const XC(func_info_type) XC(func_info_gga_x_ssb_d) = {
  XC_GGA_X_SSB_D,
  XC_EXCHANGE,
  "Swart, Sola and Bickelhaupt dispersion",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_ssb_d_init, NULL, 
  NULL, NULL, NULL
};


