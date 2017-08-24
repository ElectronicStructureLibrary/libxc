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

#define XC_HYB_GGA_XC_HSE03       427 /* the 2003 version of the screened hybrid HSE */
#define XC_HYB_GGA_XC_HSE06       428 /* the 2006 version of the screened hybrid HSE */
#define XC_HYB_GGA_XC_HJS_PBE     429 /* HJS hybrid screened exchange PBE version */
#define XC_HYB_GGA_XC_HJS_PBE_SOL 430 /* HJS hybrid screened exchange PBE_SOL version */
#define XC_HYB_GGA_XC_HJS_B88     431 /* HJS hybrid screened exchange B88 version */
#define XC_HYB_GGA_XC_HJS_B97X    432 /* HJS hybrid screened exchange B97x version */
#define XC_HYB_GGA_XC_LRC_WPBEH   465 /* Long-range corrected functional by Rorhdanz et al */
#define XC_HYB_GGA_XC_LRC_WPBE    473 /* Long-range corrected functional by Rorhdanz et al */
#define XC_HYB_GGA_XC_LC_WPBE     478 /* Long-range corrected functional by Vydrov and Scuseria */
#define XC_HYB_GGA_XC_HSE12       479 /* HSE12 by Moussa, Schultz and Chelikowsky */
#define XC_HYB_GGA_XC_HSE12S      480 /* Short-range HSE12 by Moussa, Schultz, and Chelikowsky */
#define XC_HYB_GGA_XC_HSE_SOL     481 /* HSEsol functional by Schimka, Harl, and Kresse */

static void
hyb_gga_xc_hse_init(XC(func_type) *p)
{
  /* Amount of short-range exact exchange */
  double beta;
  /* Range separation parameter */
  double omega_HF, omega_PBE;

  /* Note that there is an enormous mess in the literature concerning
     the values of omega in HSE. This is due to an error in the
     original paper that stated that they had used omega=0.15. This
     was in fact not true, and the real value used was omega^HF =
     0.15/sqrt(2) ~ 0.1061 and omega^PBE = 0.15*cbrt(2) ~ 0.1890. In
     2006 Krukau et al [JCP 125, 224106 (2006)] tried to clarify the
     situation, and called HSE03 to the above choice of parameters,
     and called HSE06 to the functional where omega^HF=omega^PBE. By
     testing several properties for atoms they reached the conclusion
     that the best value for omega=0.11.

     Of course, codes are just as messy as the papers. In espresso
     HSE06 has the value omega=0.106. VASP, on the other hand, uses
     for HSE03 the same value omega^HF = omega^PBE = 0.3 (A^-1) ~
     0.1587 and for HSE06 omega^HF = omega^PBE = 0.2 (A^-1) ~ 0.1058.

     We try to follow the original definition of the functional. The
     default omega for XC_GGA_X_WPBEH is zero, so WPBEh reduces to
     PBEh
   */
  switch(p->info->number){
  case XC_HYB_GGA_XC_HSE03:
    /* in this case one should use omega^HF = 0.15/sqrt(2) and
       omega^PBE = 0.15*CBRT(2.0)*/
    beta = 0.25;
    omega_HF = 0.15/M_SQRT2;
    omega_PBE = 0.15*CBRT(2.0);
    break;
  case XC_HYB_GGA_XC_HSE06:
    /* in this case one should use omega^HF = omega^PBE = 0.11 */
    beta = 0.25;
    omega_HF = omega_PBE = 0.11;
    break;
  case XC_HYB_GGA_XC_HSE12:
    /* omega = 0.185 with 31.3% exact short-range exchange */
    beta = 0.313;
    omega_HF = omega_PBE = 0.185;
    break;
  case XC_HYB_GGA_XC_HSE12S:
    /* omega = 0.408 with 42.5% exact short-range exchange */
    beta = 0.425;
    omega_HF = omega_PBE = 0.408;
    break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_hse\n");
    exit(1);
  }

  int   funcs_id  [3] = {XC_GGA_X_WPBEH, XC_GGA_X_WPBEH, XC_GGA_C_PBE};
  FLOAT funcs_coef[3] = {1.0, -beta, 1.0};
  
  XC(mix_init)(p, 3, funcs_id, funcs_coef);
  p->cam_omega = omega_HF;
  p->cam_beta  =  beta;
  XC(gga_x_wpbeh_set_params)(p->func_aux[1], omega_PBE);
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_hse03) = {
  XC_HYB_GGA_XC_HSE03,
  XC_EXCHANGE_CORRELATION,
  "HSE03",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Heyd2003_8207, &xc_ref_Heyd2003_8207_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hse06) = {
  XC_HYB_GGA_XC_HSE06,
  XC_EXCHANGE_CORRELATION,
  "HSE06",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Heyd2003_8207, &xc_ref_Heyd2003_8207_err, &xc_ref_Krukau2006_224106, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hse12) = {
  XC_HYB_GGA_XC_HSE12,
  XC_EXCHANGE_CORRELATION,
  "HSE12",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Moussa2012_204117, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hse12s) = {
  XC_HYB_GGA_XC_HSE12S,
  XC_EXCHANGE_CORRELATION,
  "HSE12 (short-range version)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Moussa2012_204117, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL, NULL
};

static void
hyb_gga_xc_hse_sol_init(XC(func_type) *p)
{
  int   funcs_id  [3] = {XC_GGA_X_HJS_PBE_SOL, XC_GGA_X_HJS_PBE_SOL, XC_GGA_C_PBE};
  FLOAT funcs_coef[3] = {1.0, -0.25, 1.0};
  
  XC(mix_init)(p, 3, funcs_id, funcs_coef);
  p->cam_omega = 0.11;
  p->cam_beta  = 0.25;
  XC(gga_x_hjs_set_params)(p->func_aux[0], 0.0);
  XC(gga_x_hjs_set_params)(p->func_aux[1], p->cam_omega);
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_hse_sol) = {
  XC_HYB_GGA_XC_HSE_SOL,
  XC_EXCHANGE_CORRELATION,
  "HSEsol",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schimka2011_024116, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hse_sol_init,
  NULL, NULL, NULL, NULL
};

static void
hyb_gga_xc_lc_wpbe_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_WPBEH, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0, 1.0};
  
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  
  p->cam_omega =  0.4;
  p->cam_alpha =  1.0;
  p->cam_beta  = -1.0;
  XC(gga_x_wpbeh_set_params)(p->func_aux[0], p->cam_omega);
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_lc_wpbe) = {
  XC_HYB_GGA_XC_LC_WPBE,
  XC_EXCHANGE_CORRELATION,
  "Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Vydrov2006_234109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_lc_wpbe_init,
  NULL, NULL, NULL, NULL
};


static void
hyb_gga_xc_lrc_wpbeh_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {0.8, 1.0};
  
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  
  p->cam_omega =  0.2;
  p->cam_alpha =  1.0;
  p->cam_beta  = -0.8;
  XC(gga_x_wpbeh_set_params)(p->func_aux[0], p->cam_omega);
}

static void
hyb_gga_xc_lrc_wpbe_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {1.0, 1.0};
  
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  
  p->cam_omega =  0.3;
  p->cam_alpha =  1.0;
  p->cam_beta  = -1.0;
  XC(gga_x_wpbeh_set_params)(p->func_aux[0], p->cam_omega);
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_lrc_wpbeh) = {
  XC_HYB_GGA_XC_LRC_WPBEH,
  XC_EXCHANGE_CORRELATION,
  "Long-range corrected short-range hybrid PBE (LRC-wPBEh) by Rohrdanz, Martins and Herbert",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Rohrdanz2009_054112, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_lrc_wpbeh_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_lrc_wpbe) = {
  XC_HYB_GGA_XC_LRC_WPBE,
  XC_EXCHANGE_CORRELATION,
  "Long-range corrected PBE (LRC-wPBE) by Rohrdanz, Martins and Herbert",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Rohrdanz2009_054112, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_lrc_wpbe_init,
  NULL, NULL, NULL, NULL
};

static void
hyb_gga_xc_hjs_init(XC(func_type) *p)
{
  static int   funcs_id  [3] = {-1, -1, XC_GGA_C_PBE};
  static FLOAT funcs_coef[3] = {1.0, -0.25, 1.0};

  switch(p->info->number){
  case XC_HYB_GGA_XC_HJS_PBE:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_PBE;
    break;
  case XC_HYB_GGA_XC_HJS_PBE_SOL:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_PBE_SOL;
    break;
  case XC_HYB_GGA_XC_HJS_B88:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_B88;
    break;
  case XC_HYB_GGA_XC_HJS_B97X:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_B97X;
    break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_hjs\n");
    exit(1);
  }

  XC(mix_init)(p, 3, funcs_id, funcs_coef);

  p->cam_omega = 0.11;
  p->cam_beta  = 0.25;
  XC(gga_x_hjs_set_params)(p->func_aux[0], 0.0);
  XC(gga_x_hjs_set_params)(p->func_aux[1], p->cam_omega);
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_pbe) = {
  XC_HYB_GGA_XC_HJS_PBE,
  XC_EXCHANGE_CORRELATION,
  "HJS hybrid screened exchange PBE version",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hjs_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_pbe_sol) = {
  XC_HYB_GGA_XC_HJS_PBE_SOL,
  XC_EXCHANGE_CORRELATION,
  "HJS hybrid screened exchange PBE_SOL version",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hjs_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_b88) = {
  XC_HYB_GGA_XC_HJS_B88,
  XC_EXCHANGE_CORRELATION,
  "HJS hybrid screened exchange B88 version",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hjs_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hjs_b97x) = {
  XC_HYB_GGA_XC_HJS_B97X,
  XC_EXCHANGE_CORRELATION,
  "HJS hybrid screened exchange B97x version",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_gga_xc_hjs_init,
  NULL, NULL, NULL, NULL
};
