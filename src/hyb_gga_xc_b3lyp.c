/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_HYB_GGA_XC_B3PW91        401 /* The original (ACM) hybrid of Becke    */
#define XC_HYB_GGA_XC_B3LYP         402 /* The (in)famous B3LYP                  */
#define XC_HYB_GGA_XC_B3P86         403 /* Perdew 86 hybrid similar to B3PW91    */
#define XC_HYB_GGA_XC_MPW3PW        415 /* mixture with the mPW functional       */
#define XC_HYB_GGA_XC_MPW3LYP       419 /* mixture of mPW and LYP                */
#define XC_HYB_GGA_XC_MB3LYP_RC04   437 /* B3LYP with RC04 LDA                   */
#define XC_HYB_GGA_XC_REVB3LYP      454 /* Revised B3LYP                         */
#define XC_HYB_GGA_XC_B3LYPs        459 /* B3LYP* functional                     */
#define XC_HYB_GGA_XC_B3LYP5        475 /* B3LYP with VWN functional 5 instead of RPA */
#define XC_HYB_GGA_XC_B5050LYP      572 /* Like B3LYP but more exact exchange    */
#define XC_HYB_GGA_XC_KMLYP         485 /* Kang-Musgrave hybrid                  */
#define XC_HYB_GGA_XC_APF           409 /* APF hybrid density functional         */
#define XC_HYB_GGA_XC_WC04          611 /* hybrid fitted to carbon NMR shifts    */
#define XC_HYB_GGA_XC_WP04          612 /* hybrid fitted to proton NMR shifts    */
#define XC_HYB_GGA_XC_QTP17         460 /* global hybrid for vertical ionization potentials */
#define XC_HYB_GGA_XC_B3LYP_MCM1    461 /* B3LYP reoptimized in 6-31+G(2df,p) for enthalpies of formation */
#define XC_HYB_GGA_XC_B3LYP_MCM2    462 /* B3LYP reoptimized in 6-31+G(2df,p) for enthalpies of formation */

/*************************************************************/

static const func_params_type ext_params_b3pw91[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of PW91 correlation correction"}
};

static void
set_ext_params_b3pw91(xc_func_type *p, const double *ext_params)
{
  double a0, ax, ac;

  assert(p != NULL);

  a0 = get_ext_param(p->info->ext_params, ext_params, 0);
  ax = get_ext_param(p->info->ext_params, ext_params, 1);
  ac = get_ext_param(p->info->ext_params, ext_params, 2);

  p->mix_coef[0] = 1.0 - a0 - ax;
  p->mix_coef[1] = ax;
  p->mix_coef[2] = 1.0 - ac;
  p->mix_coef[3] = ac;

  p->cam_alpha = a0;
}

void
xc_hyb_gga_xc_b3pw91_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_PW, XC_GGA_C_PW91};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3pw91 = {
  XC_HYB_GGA_XC_B3PW91,
  XC_EXCHANGE_CORRELATION,
  "The original (ACM, B3PW91) hybrid of Becke",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Becke1993_5648, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b3pw91, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3pw91_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/

static const func_params_type ext_params_b3lyp[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of LYP correlation correction"}
};

void
xc_hyb_gga_xc_b3lyp_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3lyp = {
  XC_HYB_GGA_XC_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Stephens1994_11623, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b3lyp, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3lyp_init, NULL,
  NULL, NULL, NULL
};

/*************************************************************/
void
xc_hyb_gga_xc_b3lyp5_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3lyp5 = {
  XC_HYB_GGA_XC_B3LYP5,
  XC_EXCHANGE_CORRELATION,
  "B3LYP with VWN functional 5 instead of RPA",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Stephens1994_11623, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b3lyp, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3lyp5_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/

static const func_params_type ext_params_b3p86[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of P86 correlation correction"}
};

void
xc_hyb_gga_xc_b3p86_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_P86};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3p86 = {
  XC_HYB_GGA_XC_B3P86,
  XC_EXCHANGE_CORRELATION,
  "B3P86",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_gaussianimplementation, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b3p86, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3p86_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_mpw3pw[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of mPW91 exchange correction"},
  {"_ac", 0.81, "Fraction of PW91 correlation correction"}
};

void
xc_hyb_gga_xc_mpw3pw_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_MPW91, XC_LDA_C_VWN_RPA, XC_GGA_C_PW91};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_mpw3pw = {
  XC_HYB_GGA_XC_MPW3PW,
  XC_EXCHANGE_CORRELATION,
  "MPW3PW of Adamo & Barone",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Adamo1998_664, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_mpw3pw, set_ext_params_b3pw91,
  xc_hyb_gga_xc_mpw3pw_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_mpw3lyp[] = {
  {"_a0", 0.218, "Fraction of exact exchange"},
  {"_ax", 0.709, "Fraction of mPW91 exchange correction"},
  {"_ac", 0.871, "Fraction of LYP correlation correction"}
};

void
xc_hyb_gga_xc_mpw3lyp_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_MPW91, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_mpw3lyp = {
  XC_HYB_GGA_XC_MPW3LYP,
  XC_EXCHANGE_CORRELATION,
  "MPW3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Zhao2004_6908, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_mpw3lyp, set_ext_params_b3pw91,
  xc_hyb_gga_xc_mpw3lyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_mb3lyp_rc04[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of LYP correlation correction"},
  {"_d",  0.57, "Correction factor for RC04 part"}
};

static void
set_ext_params_mb3lyp_rc04(xc_func_type *p, const double *ext_params)
{
  double a0, ax, ac, d;

  assert(p != NULL);

  a0 = get_ext_param(p->info->ext_params, ext_params, 0);
  ax = get_ext_param(p->info->ext_params, ext_params, 1);
  ac = get_ext_param(p->info->ext_params, ext_params, 2);
  d = get_ext_param(p->info->ext_params, ext_params, 3);

  p->mix_coef[0] = 1.0 - a0 - ax;
  p->mix_coef[1] = ax;
  p->mix_coef[2] = 1.0 - d*ac;
  p->mix_coef[3] = ac;

  p->cam_alpha = a0;
}

void
xc_hyb_gga_xc_mb3lyp_rc04_init(xc_func_type *p)
{
  static int    funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_RC04, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_mb3lyp_rc04 = {
  XC_HYB_GGA_XC_MB3LYP_RC04,
  XC_EXCHANGE_CORRELATION,
  "B3LYP with RC04 LDA",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Tognetti2007_381, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  4, ext_params_mb3lyp_rc04, set_ext_params_mb3lyp_rc04,
  xc_hyb_gga_xc_mb3lyp_rc04_init, NULL,
  NULL, NULL, NULL
};

/*************************************************************/

static const func_params_type ext_params_revb3lyp[] = {
  {"_a0", 0.20, "Fraction of exact exchange"},
  {"_ax", 0.67, "Fraction of B88 exchange correction"},
  {"_ac", 0.84, "Fraction of LYP correlation correction"}
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_revb3lyp = {
  XC_HYB_GGA_XC_REVB3LYP,
  XC_EXCHANGE_CORRELATION,
  "Revised B3LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Lu2013_64, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_revb3lyp, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3lyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_b3lyps[] = {
  {"_a0", 0.15, "Fraction of exact exchange"},
  {"_ax", 0.72, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of LYP correlation correction"}
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3lyps = {
  XC_HYB_GGA_XC_B3LYPs,
  XC_EXCHANGE_CORRELATION,
  "B3LYP*",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Reiher2001_48, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b3lyps, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3lyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/

static const func_params_type ext_params_b5050lyp[] = {
  {"_a0", 0.50, "Fraction of exact exchange"},
  {"_ax", 0.42, "Fraction of B88 exchange correction"},
  {"_ac", 0.81, "Fraction of LYP correlation correction"}
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b5050lyp = {
  XC_HYB_GGA_XC_B5050LYP,
  XC_EXCHANGE_CORRELATION,
  "B5050LYP",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Shao2003_4807, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  3, ext_params_b5050lyp, set_ext_params_b3pw91,
  xc_hyb_gga_xc_b3lyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_kmlyp[] = {
  {"_a0", 0.557, "Fraction of exact exchange"},
  {"_ac", 0.448, "Fraction of LYP correlation correction"},
};

static void
set_ext_params_kmlyp(xc_func_type *p, const double *ext_params)
{
  double a0, ax, ac;

  assert(p != NULL);

  a0 = get_ext_param(p->info->ext_params, ext_params, 0);
  ac = get_ext_param(p->info->ext_params, ext_params, 1);

  p->mix_coef[0] = 1.0 - a0;
  p->mix_coef[1] = 1.0 - ac;
  p->mix_coef[2] = ac;

  p->cam_alpha = a0;
}

void
xc_hyb_gga_xc_kmlyp_init(xc_func_type *p)
{
  static int    funcs_id  [3] = {XC_LDA_X, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 3, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_kmlyp = {
  XC_HYB_GGA_XC_KMLYP,
  XC_EXCHANGE_CORRELATION,
  "Kang-Musgrave hybrid",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Kang2001_11040, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  2, ext_params_kmlyp, set_ext_params_kmlyp,
  xc_hyb_gga_xc_kmlyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_qtp17[] = {
  {"_a0", 0.62, "Fraction of exact exchange"},
  {"_ac", 0.80, "Fraction of LYP correlation correction"},
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_qtp17 = {
  XC_HYB_GGA_XC_QTP17,
  XC_EXCHANGE_CORRELATION,
  "Global hybrid for vertical ionization potentials",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Jin2018_064111, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  2, ext_params_qtp17, set_ext_params_kmlyp,
  xc_hyb_gga_xc_kmlyp_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
void
xc_hyb_gga_xc_apf_init(xc_func_type *p)
{
  /* Functional is 41.1% B3PW91 and 58.9% PBE0 */
  const double fb3pw91 = 0.411;
  const double fpbe0   = 1.0 - fb3pw91;

  /* Exact exchange in B3PW91 and PBE0 */
  const double xb3pw91 = 0.20;
  const double xpbe0   = 0.25;

  int funcs_id [6] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_PW, XC_GGA_C_PW91, XC_GGA_X_PBE, XC_GGA_C_PBE};

  /* Used C standard doesn't allow initializer list, even with const
     variables */
  double funcs_coef[6];
  funcs_coef[0]=(1.0 - xb3pw91 - 0.72)*fb3pw91;
  funcs_coef[1]=0.72*fb3pw91;
  funcs_coef[2]=(1.0 - 0.81)*fb3pw91;
  funcs_coef[3]=0.81*fb3pw91;
  funcs_coef[4]=(1.0 - xpbe0)*fpbe0;
  funcs_coef[5]=fpbe0;

  xc_mix_init(p, 6, funcs_id, funcs_coef);
  p->cam_alpha = fb3pw91*xb3pw91 + fpbe0*xpbe0;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_apf = {
  XC_HYB_GGA_XC_APF,
  XC_EXCHANGE_CORRELATION,
  "APF hybrid functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Austin2012_4989, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  0, NULL, NULL,
  xc_hyb_gga_xc_apf_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
void
xc_hyb_gga_xc_wc04_init(xc_func_type *p)
{
  /* From the paper it is not clear if the LSDA is VWN or VWN_RPA.
     Due to the comparison to B3LYP I think the latter is more likely */
  const double PP[5] = {0.7400, 0.9999, 0.0001, 0.0001, 0.9999};
  static int funcs_id[4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  double funcs_coef[4];

  funcs_coef[0] = PP[2] - PP[1];
  funcs_coef[1] = PP[1];
  funcs_coef[2] = PP[4] - PP[3];
  funcs_coef[3] = PP[3];

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = PP[0];
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_wc04 = {
  XC_HYB_GGA_XC_WC04,
  XC_EXCHANGE_CORRELATION,
  "hybrid fitted to carbon NMR shifts",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Wiitala2006_1085, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  0, NULL, NULL,
  xc_hyb_gga_xc_wc04_init, NULL,
  NULL, NULL, NULL
};



/*************************************************************/
void
xc_hyb_gga_xc_wp04_init(xc_func_type *p)
{
  /* From the paper it is not clear if the LSDA is VWN or VWN_RPA.
     Due to the comparison to B3LYP I think the latter is more likely */
  const double PP[5] = {0.1189, 0.9614, 0.9999, 0.0001, 0.9999};
  static int funcs_id[4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  double funcs_coef[4];

  funcs_coef[0] = PP[2] - PP[1];
  funcs_coef[1] = PP[1];
  funcs_coef[2] = PP[4] - PP[3];
  funcs_coef[3] = PP[3];

  xc_mix_init(p, 4, funcs_id, funcs_coef);
  p->cam_alpha = PP[0];
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_wp04 = {
  XC_HYB_GGA_XC_WP04,
  XC_EXCHANGE_CORRELATION,
  "hybrid fitted to proton NMR shifts",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Wiitala2006_1085, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  0, NULL, NULL,
  xc_hyb_gga_xc_wp04_init, NULL,
  NULL, NULL, NULL
};


/*************************************************************/
static const func_params_type ext_params_b3lyp_mcm1[] = {
  {"_P1", 1.0000, "Scale factor for pure exchange"},
  {"_P2", 0.1986, "Fraction of exact exchange"},
  {"_P3", 0.6709, "Fraction of non-local exchange correction"},
  {"_P4", 0.8029, "Fraction of local exchange"},
  {"_P5", 1.1383, "Fraction of non-local correlation correction"},
  {"_P6", 0.9604, "Fraction of local correlation"}
};

static const func_params_type ext_params_b3lyp_mcm2[] = {
  {"_P1", 1.0000, "Scale factor for pure exchange"},
  {"_P2", 0.2228, "Fraction of exact exchange"},
  {"_P3", 0.7290, "Fraction of non-local exchange correction"},
  {"_P4", 0.8080, "Fraction of local exchange"},
  {"_P5", 0.9421, "Fraction of non-local correlation correction"},
  {"_P6", 0.9589, "Fraction of local correlation"}
};

static void
set_ext_params_b3lyp_mcm(xc_func_type *p, const double *ext_params)
{
  double p1, p2, p3, p4, p5, p6;

  assert(p != NULL);

  p1 = get_ext_param(p->info->ext_params, ext_params, 0);
  p2 = get_ext_param(p->info->ext_params, ext_params, 1);
  p3 = get_ext_param(p->info->ext_params, ext_params, 2);
  p4 = get_ext_param(p->info->ext_params, ext_params, 3);
  p5 = get_ext_param(p->info->ext_params, ext_params, 4);
  p6 = get_ext_param(p->info->ext_params, ext_params, 5);

  p->mix_coef[0] = p1*(p4 - p3);
  p->mix_coef[1] = p1*p3;
  p->mix_coef[2] = p6 - p5;
  p->mix_coef[3] = p5;

  p->cam_alpha = p2;
}

void
xc_hyb_gga_xc_b3lyp_mcm_init(xc_func_type *p)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_LDA_C_VWN_RPA, XC_GGA_C_LYP};
  static double funcs_coef[4] = {0.0, 0.0, 0.0, 0.0}; /* set by ext_params */

  xc_mix_init(p, 4, funcs_id, funcs_coef);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3lyp_mcm1 = {
  XC_HYB_GGA_XC_B3LYP_MCM1,
  XC_EXCHANGE_CORRELATION,
  "B3LYP-MCM1",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Caldeira2019_62, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  6, ext_params_b3lyp_mcm1, set_ext_params_b3lyp_mcm,
  xc_hyb_gga_xc_b3lyp_mcm_init, NULL,
  NULL, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b3lyp_mcm2 = {
  XC_HYB_GGA_XC_B3LYP_MCM2,
  XC_EXCHANGE_CORRELATION,
  "B3LYP-MCM2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Caldeira2019_62, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-32,
  6, ext_params_b3lyp_mcm2, set_ext_params_b3lyp_mcm,
  xc_hyb_gga_xc_b3lyp_mcm_init, NULL,
  NULL, NULL, NULL
};
