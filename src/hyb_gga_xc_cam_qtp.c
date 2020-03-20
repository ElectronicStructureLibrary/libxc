/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_CAM_QTP_00       490 /* CAM-QTP-00 */
#define XC_HYB_GGA_XC_CAM_QTP_01       482 /* CAM-QTP-01 */
#define XC_HYB_GGA_XC_CAM_QTP_02       491 /* CAM-QTP-02 */
#define XC_HYB_GGA_XC_LC_QTP           492 /* LC-QTP     */

void
xc_hyb_gga_xc_cam_qtp_init(xc_func_type *p)
{
  /* Variables that define the functionals, using Yanai's definitions */
  double alpha, beta, omega;
  double flyp;
  /* Variables translated into libxc definitions */
  double cam_alpha, cam_beta;
  /* Functional mix */
  static int funcs_id[4] = {XC_GGA_X_B88, XC_GGA_X_ITYH, XC_GGA_C_LYP, XC_LDA_C_VWN};
  double funcs_coef[4];
  int nfuncs;
  
  switch(p->info->number){
  case XC_HYB_GGA_XC_CAM_QTP_00:
    alpha = 0.54;
    beta  = 0.37;
    omega = 0.29;
    flyp = 0.8;
    break;
  case XC_HYB_GGA_XC_CAM_QTP_01:
    alpha = 0.23;
    beta  = 0.77;
    omega =  0.31;
    flyp = 0.8;
    break;
  case XC_HYB_GGA_XC_CAM_QTP_02:
    alpha = 0.28;
    beta  = 0.72;
    omega = 0.335;
    flyp = 1.0;
    break;
  case XC_HYB_GGA_XC_LC_QTP:
    alpha = 0.0;
    beta  = 1.0;
    omega = 0.475;
    flyp = 1.0;
    break;
  default:
    fprintf(stderr,"Internal error in hyb_gga_xc_cam_qtp_init.\n");
    exit(1);
  }

  /* N.B. The notation used in Yanai et al uses a different convention
     for alpha and beta. In libxc, alpha is the weight for full HF
     exchange, and beta is the weight for short-range-only exchange.

     In Yanai's convention, at short range there is a total of alpha
     exact exchange, and at long range, alpha+beta.
     
     alpha_libxc = alpha_Yanai + beta_Yanai
     beta_libxc  = - beta_Yanai
  */
  cam_alpha = alpha + beta;
  cam_beta = -beta;
    
  funcs_coef[0] = 1.0 - cam_alpha;
  funcs_coef[1] = -cam_beta;
  funcs_coef[2] = flyp;
  funcs_coef[3] = 1.0 - flyp;

  nfuncs = (flyp == 1.0) ? 3 : 4;
  xc_mix_init(p, nfuncs, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], &omega);

  p->cam_omega = omega;
  p->cam_alpha = cam_alpha;
  p->cam_beta  = cam_beta;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_00 = {
  XC_HYB_GGA_XC_CAM_QTP_00,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Verma2014_18A534, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_qtp_init,
  NULL, NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_01 = {
  XC_HYB_GGA_XC_CAM_QTP_01,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Jin2016_034107, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_qtp_init,
  NULL, NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_cam_qtp_02 = {
  XC_HYB_GGA_XC_CAM_QTP_02,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Haiduke2018_184106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_qtp_init,
  NULL, NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_qtp = {
  XC_HYB_GGA_XC_LC_QTP,
  XC_EXCHANGE_CORRELATION,
  "CAM-B3LYP retuned using ionization potentials of water",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Haiduke2018_184106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {0, NULL, NULL, NULL, NULL},
  xc_hyb_gga_xc_cam_qtp_init,
  NULL, NULL, NULL, NULL
};
