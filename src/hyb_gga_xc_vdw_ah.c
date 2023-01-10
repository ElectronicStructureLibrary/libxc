/*
 Copyright (C) 2022 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_VDW_DF_AHCX            803 /* Hybrid vdW-DF-ahcx functional to be used with vdW-DF nonlocal correlation */
#define XC_HYB_GGA_XC_VDW_DF2_AH             804 /* Hybrid vdW-DF2-ah functional to be used with vdW-DF2 nonlocal correlation */
#define XC_HYB_GGA_XC_VDW_DF2_AHBR           805 /* Hybrid vdW-DF2-ahbr functional to be used with vdW-DF2 nonlocal correlation */

static void
hyb_gga_xc_vdw_ah_init(xc_func_type *p)
{
  static int   funcs_id  [3] = {-1, -1, XC_LDA_C_PW};
  static double funcs_coef[3] = {1.0, 0.0, 1.0};

  switch(p->info->number){
  case XC_HYB_GGA_XC_VDW_DF_AHCX:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_CX13;
    break;
  case XC_HYB_GGA_XC_VDW_DF2_AH:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_RPW86;
    break;
  case XC_HYB_GGA_XC_VDW_DF2_AHBR:
    funcs_id[0] = funcs_id[1] = XC_GGA_X_HJS_B86R;
    break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_vdw_ah_init\n");
    exit(1);
  }

  xc_mix_init(p, 3, funcs_id, funcs_coef);
  xc_hyb_init_sr(p, 0.0, 0.0);
}

static void
ah_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double beta, omega;

  assert(p != NULL);

  beta  = get_ext_param(p, ext_params, 0);
  omega = get_ext_param(p, ext_params, 1);

  p->mix_coef[1] = -beta;

  p->hyb_coeff[0] = beta;
  p->hyb_omega[0] = omega;

  xc_func_set_ext_params_name(p->func_aux[0], "_omega", 0.0);
  xc_func_set_ext_params_name(p->func_aux[1], "_omega", omega);
}

#define N_PAR 2
static const double par_ahcx[N_PAR] = {0.20, 0.106};
static const double par_ah[N_PAR] = {0.25, 0.106};
static const double par_ahbr[N_PAR] = {0.25, 0.106};
static const char  *names[N_PAR] = {"_beta", "_omega"};
static const char  *desc[N_PAR]  = {
  "Mixing parameter",
  "Screening parameter"
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_vdw_df_ahcx = {
  XC_HYB_GGA_XC_VDW_DF_AHCX,
  XC_EXCHANGE_CORRELATION,
  "Hybrid vdW-DF-ahcx functional to be used with vdW-DF nonlocal correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_025902, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_ahcx, ah_set_ext_params},
  hyb_gga_xc_vdw_ah_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_vdw_df2_ah = {
  XC_HYB_GGA_XC_VDW_DF2_AH,
  XC_EXCHANGE_CORRELATION,
  "Hybrid vdW-DF2-ah functional to be used with vdW-DF2 nonlocal correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_025902, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_ah, ah_set_ext_params},
  hyb_gga_xc_vdw_ah_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_vdw_df2_ahbr = {
  XC_HYB_GGA_XC_VDW_DF2_AHBR,
  XC_EXCHANGE_CORRELATION,
  "Hybrid vdW-DF2-ahbr functional to be used with vdW-DF2 nonlocal correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_041003, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_ahbr, ah_set_ext_params},
  hyb_gga_xc_vdw_ah_init, NULL,
  NULL, NULL, NULL
};
