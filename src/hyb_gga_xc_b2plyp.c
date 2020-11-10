/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_B2PLYP    713 /* Double hybrid of Grimme */

#define B2PLYP_N_PAR 2
static const char  *b2plyp_names[B2PLYP_N_PAR]  = {"_ax", "_c"};
static const char  *b2plyp_desc[B2PLYP_N_PAR]   = {"Fock fraction", "PT2 fraction"};
static const double b2plyp_values[B2PLYP_N_PAR] = {0.53, 0.27};

static void
b2plyp_set_ext_params(xc_func_type *p, const double *ext_params)
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
hyb_gga_xc_b2plyp_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_GGA_C_LYP};
  static double funcs_coef[2] = {1.0 - 0.53, 1.0 - 0.27};

  int hyb_type[2]     = {XC_HYB_PT2, XC_HYB_FOCK};
  double hyb_omega[2] = {0.0, 0.0};
  double hyb_coeff[2] = {0.0, 0.0};

  /* Note that the value of funcs_coef[0] and hyb_coeff will be set
      by set_ext_params */
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init(p, 2, hyb_type, hyb_coeff, hyb_omega);
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b2plyp = {
  XC_HYB_GGA_XC_B2PLYP,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Grimme",
  XC_FAMILY_GGA,
  {&xc_ref_Grimme2006_034108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {B2PLYP_N_PAR, b2plyp_names, b2plyp_desc, b2plyp_values, b2plyp_set_ext_params},
  hyb_gga_xc_b2plyp_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};
