/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define  XC_HYB_GGA_XC_LC_BLYP 400  /* Long-range corrected BLYP */
#define  XC_HYB_GGA_XC_LC_BOP  636  /* Long-range corrected OP_B88 */

#define N_PAR 1
static const char *names[N_PAR] = {"_omega"};
static const char *desc[N_PAR] = {
  "Range separation parameter"
};

static const double par_lc_blyp[N_PAR] = {0.3};
static const double par_lc_bop[N_PAR] = {0.47};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  double omega;

  assert(p != NULL);
  omega  = get_ext_param(p, ext_params, 0);
  xc_func_set_ext_params_name(p->func_aux[0], "_omega", omega);

  assert(p->hyb_number_terms == 2);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = -1.0;
  p->hyb_omega[0] = omega;

  p->hyb_type[1]  = XC_HYB_FOCK;
  p->hyb_coeff[1] = 1.0;
  p->hyb_omega[1] = omega;
}


void
xc_hyb_gga_xc_lc_blyp_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_ITYH, XC_GGA_C_LYP};
  static double funcs_coef[2] = {1.0, 1.0};
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0); /* Set by parameters */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_blyp = {
  XC_HYB_GGA_XC_LC_BLYP,
  XC_EXCHANGE_CORRELATION,
  "LC version of BLYP",
  XC_FAMILY_GGA,
  {&xc_ref_Anderson2017_1656, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {N_PAR, names, desc, par_lc_blyp, set_ext_params},
  xc_hyb_gga_xc_lc_blyp_init, NULL,
  NULL, NULL, NULL
};

void
xc_hyb_gga_xc_lc_bop_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_ITYH, XC_GGA_C_OP_B88};
  static double funcs_coef[2] = {1.0, 1.0};
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0); /* set by parameters */
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_bop = {
  XC_HYB_GGA_XC_LC_BOP,
  XC_EXCHANGE_CORRELATION,
  "LC version of OP_B88",
  XC_FAMILY_GGA,
  {&xc_ref_Song2007_154105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-14,
  {N_PAR, names, desc, par_lc_bop, set_ext_params},
  xc_hyb_gga_xc_lc_bop_init, NULL,
  NULL, NULL, NULL
};
