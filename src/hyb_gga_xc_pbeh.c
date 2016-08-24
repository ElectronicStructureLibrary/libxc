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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_HYB_GGA_XC_PBEH      406 /* aka PBE0 or PBE1PBE */
#define XC_HYB_GGA_XC_PBE0_13   456 /* PBE0-1/3            */
#define XC_HYB_GGA_XC_HPBEINT   472 /* hPBEint             */
#define XC_HYB_GGA_XC_PBE_MOL0  273 /* PBEmol0             */
#define XC_HYB_GGA_XC_PBE_SOL0  274 /* PBEsol0             */
#define XC_HYB_GGA_XC_PBEB0     275 /* PBEbeta0            */
#define XC_HYB_GGA_XC_PBE_MOLB0 276 /* PBEmolbeta0         */

static void
hyb_gga_xc_pbeh_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

void 
XC(hyb_gga_xc_pbeh_set_params)(XC(func_type) *p, FLOAT alpha)
{
  assert(alpha>=0 && alpha<=1.0);

  p->cam_alpha   = alpha;
  p->mix_coef[0] = 1.0 - alpha;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbeh) = {
  XC_HYB_GGA_XC_PBEH,
  XC_EXCHANGE_CORRELATION,
  "PBEH (PBE0)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Adamo1999_6158, &xc_ref_Ernzerhof1999_5029, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbeh_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

static void
hyb_gga_xc_pbe0_13_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0 - 1.0/3.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 1.0/3.0;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbe0_13) = {
  XC_HYB_GGA_XC_PBE0_13,
  XC_EXCHANGE_CORRELATION,
  "PBE0-1/3",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Cortona2012_086101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbe0_13_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_gga_xc_hpbeint_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBEINT, XC_GGA_C_PBEINT};
  static FLOAT funcs_coef[2] = {1.0 - 1.0/6.0, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 1.0/6.0;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_hpbeint) = {
  XC_HYB_GGA_XC_HPBEINT,
  XC_EXCHANGE_CORRELATION,
  "hPBEint",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Fabiano2013_673, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_DEVELOPMENT,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hpbeint_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_gga_xc_pbemol0_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE_MOL, XC_GGA_C_PBE_MOL};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbe_mol0) = {
  XC_HYB_GGA_XC_PBE_MOL0,
  XC_EXCHANGE_CORRELATION,
  "PBEmol0",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbemol0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_gga_xc_pbesol0_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE_SOL, XC_GGA_C_PBE_SOL};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbe_sol0) = {
  XC_HYB_GGA_XC_PBE_SOL0,
  XC_EXCHANGE_CORRELATION,
  "PBEsol0",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbesol0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_gga_xc_pbeb0_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  /* 0.050044 ~ 3/4 beta_PBE */
  XC(gga_c_pbe_set_params)(p->func_aux[1], 0.050044);
  p->cam_alpha = 0.25;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbeb0) = {
  XC_HYB_GGA_XC_PBEB0,
  XC_EXCHANGE_CORRELATION,
  "PBEbeta0",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbeb0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};


static void
hyb_gga_xc_pbemolb0_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE_MOL, XC_GGA_C_PBE_MOL};
  static FLOAT funcs_coef[2] = {1.0 - 0.25, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  /* 0.06288 ~ 3/4 beta_PBEmol */
  XC(gga_c_pbe_set_params)(p->func_aux[1], 0.06288);
  p->cam_alpha = 0.25;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_pbe_molb0) = {
  XC_HYB_GGA_XC_PBE_MOLB0,
  XC_EXCHANGE_CORRELATION,
  "PBEmolbeta0",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_pbemolb0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
