/*
 Copyright (C) 2015 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_XC_VV10         255 /* Vydrov and Van Voorhis */
#define XC_HYB_GGA_XC_LC_VV10  469 /* Vydrov and Van Voorhis */

#define N_PAR 2
static const char  *names[N_PAR] = {"_b", "_C" };
static const char  *desc[N_PAR]  = {
  "VV10 b parameter",
  "VV10 C parameter"
};

static const double par_vv10[N_PAR] = {5.9, 0.0093};

static void
gga_xc_vv10_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_RPW86, XC_GGA_C_PBE};
  static double funcs_coef[2] = {1.0, 1.0};
  
  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_vdw_vv10(p, 0.0, 0.0);
}

static void
vv10_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);

  /* Non-local correlation part */
  p->hyb_params[0][0] = get_ext_param(p, ext_params, 0);
  p->hyb_params[0][1] = get_ext_param(p, ext_params, 1);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_vv10 = {
  XC_GGA_XC_VV10,
  XC_EXCHANGE_CORRELATION,
  "Vydrov and Van Voorhis",
  XC_FAMILY_GGA,
  {&xc_ref_Vydrov2010_244103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, par_vv10, vv10_set_ext_params},
  gga_xc_vv10_init,
  NULL, NULL, NULL, NULL
};

#define LC_N_PAR 5
static const char  *lc_names[LC_N_PAR] = {"_alpha", "_beta", "_omega", "_b", "_C" };
static const char  *lc_desc[LC_N_PAR]  = {
  "Fraction of Hartree-Fock exchange",
  "Fraction of short-range exact exchange",
  "Range separation constant",
  "VV10 b parameter",
  "VV10 C parameter"
};

static const double par_lc_vv10[LC_N_PAR] = {1.00, -1.00, 0.45, 6.3, 0.0089};

static void
hyb_gga_xc_lc_vv10_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_HJS_PBE, XC_GGA_C_PBE};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);

  p->hyb_number_terms = 3; /* we add a vv10 term */
  p->hyb_type[2] = XC_HYB_VDW_VV10;
}

static void
lc_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);

  p->hyb_params[0][0] = get_ext_param(p, ext_params, 0); /* alpha */
  p->hyb_params[1][0] = get_ext_param(p, ext_params, 1); /* beta  */
  p->hyb_params[1][1] = get_ext_param(p, ext_params, 2); /* omega */
  p->hyb_params[2][0] = get_ext_param(p, ext_params, 3); /* b */
  p->hyb_params[2][1] = get_ext_param(p, ext_params, 4); /* C */

  /* DFT part */
  p->mix_coef[0] = -p->hyb_params[1][0];
  xc_func_set_ext_params_name(p->func_aux[0], "_omega", p->hyb_params[1][1]);
}


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_lc_vv10 = {
  XC_HYB_GGA_XC_LC_VV10,
  XC_EXCHANGE_CORRELATION,
  "Vydrov and Van Voorhis",
  XC_FAMILY_GGA,
  {&xc_ref_Vydrov2010_244103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {LC_N_PAR, lc_names, lc_desc, par_lc_vv10, lc_set_ext_params},
  hyb_gga_xc_lc_vv10_init,
  NULL, NULL, NULL, NULL
};
