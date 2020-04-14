/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_SRC1_BLYP       714 /* Hybrid with two range separations (form 1) */
#define XC_HYB_GGA_XC_SRC2_BLYP       715 /* Hybrid with two range separations (form 2) */

#define SRC1_N_PAR 6
static const char  *src1_names[SRC1_N_PAR] = {"_C_SR", "_omega_SR", "_C_LR", "_omega_LR", "_C_LYP", "_C_VWN"};
static const char  *src1_desc[SRC1_N_PAR]  = {
  "Short-range mixing parameter",
  "Short-range screening parameter",
  "Long-range mixing parameter",
  "Long-range screening parameter",
  "coefficient of LYP",
  "coefficient of VWN",
};
static const double src1_values[SRC1_N_PAR]  = {0.50, 0.56, 0.17, 2.45, 0.81, 0.19};
static const double src2_values[SRC1_N_PAR]  = {0.55, 0.69, 0.08, 1.02, 0.81, 0.19};

static void
hyb_gga_xc_src1_init(xc_func_type *p)
{
  int    funcs_id  [5] = {XC_GGA_X_B88, XC_GGA_X_ITYH, XC_GGA_X_ITYH, XC_GGA_C_LYP, XC_LDA_C_VWN};
  double funcs_coef[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  int    hyb_type[]  = {XC_HYB_ERF_SR, XC_HYB_ERF_SR, XC_HYB_FOCK};
  double hyb_coeff[] = {0.0, 0.0, 0.0};
  double hyb_omega[] = {0.0, 0.0, 0.0};
  
  /* Note that the value of funcs_coef will be set by set_ext_params */
  xc_mix_init(p, 5, funcs_id, funcs_coef);
  xc_hyb_init(p, 3, hyb_type, hyb_coeff, hyb_omega);
}

static void
src1_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double C_SR, C_LR;

  assert(p != NULL);

  C_SR            = get_ext_param(p, ext_params, 0);
  p->hyb_omega[0] = get_ext_param(p, ext_params, 1); /* omega SR */
  C_LR            = get_ext_param(p, ext_params, 2);
  p->hyb_omega[1] = get_ext_param(p, ext_params, 3); /* omega_LR */
  p->mix_coef[3]  = get_ext_param(p, ext_params, 4); /* LYP mixing */
  p->mix_coef[4]  = get_ext_param(p, ext_params, 5); /* VWN mixing */

  switch(p->info->number){
  case XC_HYB_GGA_XC_SRC1_BLYP:
    p->mix_coef[0] = 1.0 - C_LR;
    p->mix_coef[1] =-C_SR;
    p->mix_coef[2] = C_LR;
    break;
  case XC_HYB_GGA_XC_SRC2_BLYP:
    p->mix_coef[0] = 1.0 - C_LR;
    p->mix_coef[1] = 1.0 - C_SR;
    p->mix_coef[2] = C_LR - 1.0;
    break;
  }
  p->hyb_coeff[0] = C_SR;
  p->hyb_coeff[1] =-C_LR;
  p->hyb_coeff[2] = C_LR; /* Normal Fock */
  
  xc_func_set_ext_params_name(p->func_aux[1], "_omega", p->hyb_omega[0]); /* mu_SR */
  xc_func_set_ext_params_name(p->func_aux[2], "_omega", p->hyb_omega[1]); /* mu_LR */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_src1_blyp = {
  XC_HYB_GGA_XC_SRC1_BLYP,
  XC_EXCHANGE_CORRELATION,
  "Hybrid with two range separations (form 1)",
  XC_FAMILY_GGA,
  {&xc_ref_Besley2009_10350, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {SRC1_N_PAR, src1_names, src1_desc, src1_values, src1_set_ext_params},
  hyb_gga_xc_src1_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_src2_blyp = {
  XC_HYB_GGA_XC_SRC2_BLYP,
  XC_EXCHANGE_CORRELATION,
  "Hybrid with two range separations (form 2)",
  XC_FAMILY_GGA,
  {&xc_ref_Besley2009_10350, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  {SRC1_N_PAR, src1_names, src1_desc, src2_values, src1_set_ext_params},
  hyb_gga_xc_src1_init, NULL,
  NULL, NULL, NULL
};

