/*
 Copyright (C) 2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_PBE0_DH   725 /* Double hybrid by Bremond and Adamo */
#define XC_HYB_GGA_XC_PBE0_2    726 /* Double hybrid by Chai and Mao */
#define XC_HYB_GGA_XC_PBE_QIDH  727 /* Double hybrid by Bremond et al */
#define XC_HYB_GGA_XC_LS1DH_PBE 728 /* Double hybrid by Toulouse et al */

#define N_PAR 2
static const char  *names[N_PAR]  = {"_ax", "_c"};
static const char  *desc[N_PAR]   = {"Fock fraction", "PT2 fraction"};
static const double pbe0_dh_values[N_PAR]   = {0.50, 0.125};
static const double pbe0_2_values[N_PAR]    = {0.7937005259840997373758528196361541301957, 0.50}; /* {cbrt(0.50), 0.50}; */
static const double pbe_qidh_values[N_PAR]  = {0.6933612743506347048433522747859617954459, 1.0/3.0}; /* {1.0/cbrt(3.0), 1.0/3.0}; */
static const double ls1dh_pbe_values[N_PAR] = {0.75, 0.75*0.75*0.75}; /* {0.75, pow(0.75,3)}; */

static void
pbe0_dh_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double ax, c;

  assert(p != NULL);
  ax = get_ext_param(p, ext_params, 0);
  c  = get_ext_param(p, ext_params, 1);

  p->mix_coef[0] = 1.0 - ax;
  p->mix_coef[1] = 1.0 - c;

  p->hyb_coeff[0] = ax;
  p->hyb_coeff[1] = c;
}

static void
hyb_gga_xc_pbe0_dh_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE, XC_GGA_C_PBE};
  static double funcs_coef[2] = {0.0, 0.0};

  int hyb_type[2]     = {XC_HYB_FOCK, XC_HYB_PT2};
  double hyb_omega[2] = {0.0, 0.0};
  double hyb_coeff[2] = {0.0, 0.0};

  /* Note that the values of funcs_coef and hyb_coeff will be set
      by set_ext_params */
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init(p, 2, hyb_type, hyb_coeff, hyb_omega);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_pbe0_dh = {
  XC_HYB_GGA_XC_PBE0_DH,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Bremond and Adamo",
  XC_FAMILY_GGA,
  {&xc_ref_Bremond2011_024106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, pbe0_dh_values, pbe0_dh_set_ext_params},
  hyb_gga_xc_pbe0_dh_init, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_pbe0_2 = {
  XC_HYB_GGA_XC_PBE0_2,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Chai and Mao",
  XC_FAMILY_GGA,
  {&xc_ref_Chai2012_121, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, pbe0_2_values, pbe0_dh_set_ext_params},
  hyb_gga_xc_pbe0_dh_init, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_pbe_qidh = {
  XC_HYB_GGA_XC_PBE_QIDH,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Bremond et al",
  XC_FAMILY_GGA,
  {&xc_ref_Bremond2014_031101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, pbe_qidh_values, pbe0_dh_set_ext_params},
  hyb_gga_xc_pbe0_dh_init, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_ls1dh_pbe = {
  XC_HYB_GGA_XC_LS1DH_PBE,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Toulouse et al",
  XC_FAMILY_GGA,
  {&xc_ref_Toulouse2011_101102, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, ls1dh_pbe_values, pbe0_dh_set_ext_params},
  hyb_gga_xc_pbe0_dh_init, NULL, NULL
};
