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

#define XC_GGA_XC_XLYP       166  /* XLYP functional */
#define XC_GGA_XC_PBE1W      173  /* Functionals fitted for water */
#define XC_GGA_XC_MPWLYP1W   174  /* Functionals fitted for water */
#define XC_GGA_XC_PBELYP1W   175  /* Functionals fitted for water */

static void
gga_xc_xlyp_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_PW91, XC_GGA_C_LYP};
  static double funcs_coef[4] = {1.0 - 0.722 - 0.347, 0.722, 0.347, 1.0};

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

const xc_func_info_type xc_func_info_gga_xc_xlyp = {
  XC_GGA_XC_XLYP,
  XC_EXCHANGE_CORRELATION,
  "XLYP",
  XC_FAMILY_GGA,
  {&xc_ref_Xu2004_2673, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_xc_xlyp_init, 
  NULL, NULL, NULL, NULL
};


static void
gga_xc_pbe1w_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_C_VWN, XC_GGA_X_PBE, XC_GGA_C_PBE};
  static double funcs_coef[3] = {1.0 - 74.0/100.0, 1.0, 74.0/100.0};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

const xc_func_info_type xc_func_info_gga_xc_pbe1w = {
  XC_GGA_XC_PBE1W,
  XC_EXCHANGE_CORRELATION,
  "PBE1W",
  XC_FAMILY_GGA,
  {&xc_ref_Dahlke2005_15677, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_xc_pbe1w_init, 
  NULL, NULL, NULL, NULL
};


static void
gga_xc_mpwlyp1w_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_C_VWN, XC_GGA_X_MPW91, XC_GGA_C_LYP};
  static double funcs_coef[3] = {1.0 - 88.0/100.0, 1.0, 88.0/100.0};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

const xc_func_info_type xc_func_info_gga_xc_mpwlyp1w = {
  XC_GGA_XC_MPWLYP1W,
  XC_EXCHANGE_CORRELATION,
  "mPWLYP1w",
  XC_FAMILY_GGA,
  {&xc_ref_Dahlke2005_15677, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_xc_mpwlyp1w_init, 
  NULL, NULL, NULL, NULL
};


static void
gga_xc_pbelyp1w_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_LDA_C_VWN, XC_GGA_X_PBE, XC_GGA_C_LYP};
  static double funcs_coef[3] = {1.0 - 74.0/100.0, 1.0, 74.0/100.0};

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

const xc_func_info_type xc_func_info_gga_xc_pbelyp1w = {
  XC_GGA_XC_PBELYP1W,
  XC_EXCHANGE_CORRELATION,
  "PBELYP1W",
  XC_FAMILY_GGA,
  {&xc_ref_Dahlke2005_15677, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_xc_pbelyp1w_init, 
  NULL, NULL, NULL, NULL
};

