/*
 Copyright (C) 2021 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_XC_VDW_DF1   641 /* original vdw_df functional of Dion et al */

#define VDW_DF1_N_PAR 1
static const char  *vdw_df1_names[VDW_DF1_N_PAR]  = {"_Z_ab"};
static const char  *vdw_df1_desc[VDW_DF1_N_PAR]   = {"parameters that enter the screened exchange"};

static const double vdw_df1_values[VDW_DF1_N_PAR] = {VDW_DF1_ZAB};

void
vdw_df1_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PBE_R, XC_LDA_C_PW};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_vdw_df(p, 0.0, 0.0); /* set by set_ext_params */
}


static void
vdw_df1_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);

  p->hyb_params[0].df.delta = 1.0;
  p->hyb_params[0].df.Zab   = get_ext_param(p, ext_params, 0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_vdw_df1 = {
  XC_GGA_XC_VDW_DF1,
  XC_EXCHANGE_CORRELATION,
  "original vdw_df functional of Dion et al",
  XC_FAMILY_GGA,
  {&xc_ref_Dion2004_246401, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {VDW_DF1_N_PAR, vdw_df1_names, vdw_df1_desc, vdw_df1_values, vdw_df1_set_ext_params},
  vdw_df1_init, NULL,
  NULL, NULL, NULL
};
