/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_CAM_B3LYP        433 /* CAM version of B3LYP */
#define XC_HYB_GGA_XC_CAMH_B3LYP       614 /* CAM version of B3LYP tuned for tddft */
#define XC_HYB_GGA_XC_TUNED_CAM_B3LYP  434 /* CAM version of B3LYP tuned for excitations */
#define XC_HYB_GGA_XC_RCAM_B3LYP       610 /* Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction */
#define XC_HYB_GGA_XC_CAM_PBEH         681 /* CAM version of PBEH */

void
xc_hyb_gga_xc_cam_b3lyp_init(xc_func_type *p)
{
  double ac = 0.81;
  static int   funcs_id  [4] = {XC_GGA_X_B88, XC_GGA_X_ITYH, XC_LDA_C_VWN, XC_GGA_C_LYP};
  double funcs_coef[4];

  /* Need temp variables since cam_ parameters are initialized in mix_init */
  static double omega, alpha, beta;

  switch(p->info->number){
  case XC_HYB_GGA_XC_CAM_B3LYP:
    /* N.B. The notation used in Yanai et al uses a different
       convention for alpha and beta.  In libxc, alpha is the weight
       for HF exchange, which in Yanai et al is alpha+beta, so:

       alpha_libxc = alpha_Yanai + beta_Yanai
       beta_libxc  = - beta_Yanai
     */
    omega = 0.33;
    alpha = 0.65;
    beta  =-0.46;
    break;
  case XC_HYB_GGA_XC_CAMH_B3LYP:
    /* The same note applies here. */
    omega = 0.33;
    alpha = 0.50;
    beta  = -0.31;
    break;
  case XC_HYB_GGA_XC_TUNED_CAM_B3LYP:
    /* The same note applies here. */
    omega = 0.150;
    alpha = 1.0000;
    beta  =-0.9201;
    break;
  default:
    fprintf(stderr,"Internal error in hyb_gga_xc_cam_b3lyp_init.\n");
    exit(1);
  }

  funcs_coef[0] = 1.0 - alpha;
  funcs_coef[1] = -beta;
  funcs_coef[2] = 1.0 - ac;
  funcs_coef[3] = ac;

  xc_mix_init(p, 4, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], &omega);

  xc_hyb_init_cam(p, omega, alpha, beta);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_b3lyp = {
  XC_HYB_GGA_XC_CAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP",
  XC_FAMILY_GGA,
  {&xc_ref_Yanai2004_51, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-9,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_camh_b3lyp = {
  XC_HYB_GGA_XC_CAMH_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP, tuned for TDDFT",
  XC_FAMILY_GGA,
  {&xc_ref_Shao2020_587, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-9,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_tuned_cam_b3lyp = {
  XC_HYB_GGA_XC_TUNED_CAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "CAM version of B3LYP, tuned for excitations and properties",
  XC_FAMILY_GGA,
  {&xc_ref_Okuno2012_29, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-9,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_b3lyp_init, NULL,
  NULL, NULL, NULL
};


void
xc_hyb_gga_xc_rcam_b3lyp_init(xc_func_type *p)
{
  static int funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_ITYH, XC_GGA_C_LYP};
  static double funcs_coef[4];

  /* Temp for cam_ parameters */
  static double omega, alpha, beta;
  /* Temp parameters for functional */
  static double a, b, cb88;

  a = 0.18352;
  b = 0.94979;
  omega = 0.33;
  cb88  = 0.95238;

  funcs_coef[0] = 1.0 - a - cb88;
  funcs_coef[1] = cb88 - b;
  funcs_coef[2] = b;
  funcs_coef[3] = 1.0;

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  xc_func_set_ext_params(p->func_aux[2], &omega);

  /* Libxc hybrid parameters */
  alpha = a + b;
  beta  =-b;

  xc_hyb_init_cam(p, omega, alpha, beta);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_rcam_b3lyp = {
  XC_HYB_GGA_XC_RCAM_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2007_191109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-9,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_rcam_b3lyp_init, NULL,
  NULL, NULL, NULL
};

#define CAM_N_PAR 4
static const char  *cam_names[CAM_N_PAR]  = {"_alpha", "_beta", "_omega_HF", "_omega_PBE"};
static const char  *cam_desc[CAM_N_PAR]   = {
  "Mixing parameter",
  "Mixing parameter in the SR",
  "Screening parameter for HF",
  "Screening parameter for PBE"
};
static const double cam_values[CAM_N_PAR] = {0.2, 0.8, 0.7, 0.7};

static void
cam_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double alpha, beta, omega_HF, omega_PBE;

  assert(p != NULL);

  alpha     = get_ext_param(p, ext_params, 0);
  beta      = get_ext_param(p, ext_params, 1);
  omega_HF  = get_ext_param(p, ext_params, 2);
  omega_PBE = get_ext_param(p, ext_params, 3);

  p->mix_coef[0] = 1.0 - alpha;
  p->mix_coef[1] = -beta;

  xc_hyb_init_cam(p, omega_HF, alpha, beta);

  xc_func_set_ext_params_name(p->func_aux[1], "_omega", omega_PBE);
}

static void
hyb_gga_xc_cam_pbeh_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {XC_GGA_X_PBE, XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static double funcs_coef[3] = {0.0, 0.0, 1.0}; /* the first two get set by ext_params */

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_pbeh = {
  XC_HYB_GGA_XC_CAM_PBEH,
  XC_EXCHANGE_CORRELATION,
  "CAM hybrid screened exchange PBE version",
  XC_FAMILY_GGA,
  {&xc_ref_Chen2018_073803, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15, 
  {CAM_N_PAR, cam_names, cam_desc, cam_values, cam_set_ext_params},
  hyb_gga_xc_cam_pbeh_init,
  NULL, NULL, NULL, NULL
};
