/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_B2PLYP     713 /* Double hybrid of Grimme */
#define XC_HYB_GGA_XC_B2GPPLYP   721 /* Double hybrid of Karton et al */
#define XC_HYB_GGA_XC_WB2PLYP    722 /* Double hybrid of Casanova-Paez, Dardis and Goerigk */
#define XC_HYB_GGA_XC_WB2GPPLYP  723 /* Double hybrid of Casanova-Paez, Dardis and Goerigk */

#define N_PAR 2
static const char  *names[N_PAR]  = {"_ax", "_c"};
static const char  *desc[N_PAR]   = {"Fock fraction", "PT2 fraction"};
static const double b2plyp_values[N_PAR]   = {0.53, 0.27};
static const double b2gpplyp_values[N_PAR] = {0.65, 0.36};

static void
hyb_gga_xc_b2plyp_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_GGA_C_LYP};
  static double funcs_coef[2] = {0.0, 0.0};

  /* Note that the values of funcs_coef and hyb_coeff will be set
      by set_ext_params */
  xc_mix_init(p, 2, funcs_id, funcs_coef);

  p->hyb_number_terms = 2;
  p->hyb_type[0] = XC_HYB_FOCK;
  p->hyb_type[1] = XC_HYB_PT2;
}

static void
b2plyp_set_ext_params(xc_func_type *p, const double *ext_params)
{
  assert(p != NULL);
  p->hyb_params[0][0] = get_ext_param(p, ext_params, 0);
  p->hyb_params[1][0] = get_ext_param(p, ext_params, 1);

  p->mix_coef[0] = 1.0 - p->hyb_params[0][0];
  p->mix_coef[1] = 1.0 - p->hyb_params[1][0];
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
  {N_PAR, names, desc, b2plyp_values, b2plyp_set_ext_params},
  hyb_gga_xc_b2plyp_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b2gpplyp = {
  XC_HYB_GGA_XC_B2GPPLYP,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Karton et al",
  XC_FAMILY_GGA,
  {&xc_ref_Karton2008_12868, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PAR, names, desc, b2gpplyp_values, b2plyp_set_ext_params},
  hyb_gga_xc_b2plyp_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};


#define N_PARW 3
static const char  *wnames[N_PARW]  = {"_ax", "_c", "_omega"};
static const char  *wdesc[N_PARW]   = {"Fock fraction", "PT2 fraction", "Range separation parameter"};
static const double wb2plyp_values[N_PARW]   = {0.53, 0.27, 0.30};
static const double wb2gpplyp_values[N_PARW] = {0.65, 0.36, 0.27};

static void
hyb_gga_xc_wb2plyp_init(xc_func_type *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_ITYH, XC_GGA_C_LYP};
  static double funcs_coef[2] = {0.0, 0.0};

  /* Note that the values of funcs_coef and hyb_coeff will be set
      by set_ext_params */
  xc_mix_init(p, 2, funcs_id, funcs_coef);

  p->hyb_number_terms = 3;
  p->hyb_type[0] = XC_HYB_FOCK;
  p->hyb_type[1] = XC_HYB_PT2;
  p->hyb_type[2] = XC_HYB_ERF_SR;
}

static void
wb2plyp_set_ext_params(xc_func_type *p, const double *ext_params)
{
  double ax, c, omega;

  assert(p != NULL);
  ax    = get_ext_param(p, ext_params, 0);
  c     = get_ext_param(p, ext_params, 1);
  omega = get_ext_param(p, ext_params, 2);

  /* Range-separation parameter */
  xc_func_set_ext_params_name(p->func_aux[0], "_omega", omega);
  p->hyb_params[2][1] = omega;

  /* Long-range coefficient is always 1 */
  p->hyb_params[0][0] = 1.0;
  p->hyb_params[1][0] = c;
  p->hyb_params[2][0] = ax - 1.0;

  /* Mixing coefficients */
  p->mix_coef[0] = 1.0 - ax;
  p->mix_coef[1] = 1.0 - c;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_wb2plyp = {
  XC_HYB_GGA_XC_WB2PLYP,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Casanova-Paez, Dardis and Goerigk",
  XC_FAMILY_GGA,
  {&xc_ref_CasanovaPaez2019_4735, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PARW, wnames, wdesc, wb2plyp_values, wb2plyp_set_ext_params},
  hyb_gga_xc_wb2plyp_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_wb2gpplyp = {
  XC_HYB_GGA_XC_WB2GPPLYP,
  XC_EXCHANGE_CORRELATION,
  "Double hybrid of Casanova-Paez, Dardis and Goerigk",
  XC_FAMILY_GGA,
  {&xc_ref_CasanovaPaez2019_4735, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-15,
  {N_PARW, wnames, wdesc, wb2gpplyp_values, wb2plyp_set_ext_params},
  hyb_gga_xc_wb2plyp_init, NULL,
  NULL, NULL, NULL /* this is taken care of by the generic routine */
};
