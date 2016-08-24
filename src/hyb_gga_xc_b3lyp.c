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

#define XC_HYB_GGA_XC_B3PW91        401 /* The original (ACM) hybrid of Becke    */
#define XC_HYB_GGA_XC_B3LYP         402 /* The (in)famous B3LYP                  */
#define XC_HYB_GGA_XC_B3P86         403 /* Perdew 86 hybrid similar to B3PW91    */
#define XC_HYB_GGA_XC_MPW3PW        415 /* mixture with the mPW functional       */
#define XC_HYB_GGA_XC_MPW3LYP       419 /* mixture of mPW and LYP                */
#define XC_HYB_GGA_XC_MB3LYP_RC04   437 /* B3LYP with RC04 LDA                   */
#define XC_HYB_GGA_XC_REVB3LYP      454 /* Revised B3LYP                         */
#define XC_HYB_GGA_XC_B3LYPs        459 /* B3LYP* functional                     */
#define XC_HYB_GGA_XC_B3LYP5        475 /* B3LYP with VWN functional 5 instead of RPA */

/*************************************************************/
void
XC(hyb_gga_xc_b3pw91_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_PW, XC_GGA_C_PW91};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3pw91) = {
  XC_HYB_GGA_XC_B3PW91,
  XC_EXCHANGE_CORRELATION,
  "The original (ACM, B3PW91) hybrid of Becke",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Becke1993_5648, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_b3pw91_init),
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_b3lyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyp) = {
  XC_HYB_GGA_XC_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Stephens1994_11623, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_b3lyp_init),
  NULL, NULL, NULL, NULL
};

/*************************************************************/
void
XC(hyb_gga_xc_b3lyp5_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyp5) = {
  XC_HYB_GGA_XC_B3LYP5,
  XC_EXCHANGE_CORRELATION,
  "B3LYP with VWN functional 5 instead of RPA",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Stephens1994_11623, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_b3lyp5_init),
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_b3p86_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_P86};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3p86) = {
  XC_HYB_GGA_XC_B3P86,
  XC_EXCHANGE_CORRELATION,
  "B3P86",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_gaussianimplementation, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_b3p86_init),
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_mpw3pw_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_MPW91, XC_LDA_C_VWN_RPA, XC_GGA_C_PW91};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3pw) = {
  XC_HYB_GGA_XC_MPW3PW,
  XC_EXCHANGE_CORRELATION,
  "MPW3PW of Adamo & Barone",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Adamo1998_664, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_mpw3pw_init), 
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_mpw3lyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_MPW91, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.218 - 0.709, 0.709, 1.0 - 0.871, 0.871};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.218;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3lyp) = {
  XC_HYB_GGA_XC_MPW3LYP,
  XC_EXCHANGE_CORRELATION,
  "MPW3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Zhao2004_6908, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_mpw3lyp_init), 
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_mb3lyp_rc04_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_RC04, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.72, 0.72, 1.0 - 0.57*0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_mb3lyp_rc04) = {
  XC_HYB_GGA_XC_MB3LYP_RC04,
  XC_EXCHANGE_CORRELATION,
  "B3LYP with RC04 LDA",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Tognetti2007_381, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_mb3lyp_rc04_init),
  NULL, NULL, NULL, NULL
};

/*************************************************************/
void
XC(hyb_gga_xc_revb3lyp_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.20 - 0.67, 0.67, 1.0 - 0.84, 0.84};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_revb3lyp) = {
  XC_HYB_GGA_XC_REVB3LYP,
  XC_EXCHANGE_CORRELATION,
  "Revised B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Lu2013_64, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_revb3lyp_init),
  NULL, NULL, NULL, NULL
};


/*************************************************************/
void
XC(hyb_gga_xc_b3lyps_init)(XC(func_type) *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.0 - 0.15 - 0.72, 0.72, 1.0 - 0.81, 0.81};

  XC(mix_init)(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = 0.15;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyps) = {
  XC_HYB_GGA_XC_B3LYPs,
  XC_EXCHANGE_CORRELATION,
  "B3LYP*",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Reiher2001_48, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  XC(hyb_gga_xc_b3lyps_init),
  NULL, NULL, NULL, NULL
};
