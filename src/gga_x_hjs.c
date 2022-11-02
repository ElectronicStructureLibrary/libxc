/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_HJS_PBE     525 /* HJS screened exchange PBE version */
#define XC_GGA_X_HJS_PBE_SOL 526 /* HJS screened exchange PBE_SOL version */
#define XC_GGA_X_HJS_B88     527 /* HJS screened exchange B88 version */
#define XC_GGA_X_HJS_B97X    528 /* HJS screened exchange B97x version */
#define XC_GGA_X_HJS_CX13    800 /* HJS screened exchange for cx13 and thus vdW-DF-ahcx */
#define XC_GGA_X_HJS_PW86    801 /* HJS screened exchange for rPW86 and thus vdW-DF2-ah */
#define XC_GGA_X_HJS_B86R    802 /* HJS screened exchange for b86r and thus vdW-DF2-ahbr */


typedef struct{
  double a[6], b[9]; /* pointers to the a and b parameters */
} gga_x_hjs_params;

#define N_PARS 16
static const char  *names[N_PARS]  = {"_a0", "_a1", "_a2", "_a3", "_a4", "_a5", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5", "_b6", "_b7", "_b8", "_omega"};
static const char  *desc[N_PARS]   = {"a0", "a1", "a2", "a3", "a4", "a5", "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "omega"};

static const double pars_PBE[N_PARS] =
  {0.0159941, 0.0852995, -0.160368, 0.152645, -0.0971263, 0.0422061,
   5.33319, -12.4780, 11.0988, -5.11013, 1.71468, -0.610380, 0.307555, -0.0770547, 0.0334840, 0.11};

static const double pars_PBE_sol[N_PARS] =
   {0.0047333, 0.0403304, -0.0574615, 0.0435395, -0.0216251, 0.0063721,
    8.52056, -13.9885, 9.28583, -3.27287, 0.843499, -0.235543, 0.0847074, -0.0171561, 0.0050552, 0.11};

static const double pars_B88[N_PARS] =
  {0.00968615, -0.0242498, 0.0259009, -0.0136606, 0.00309606, -7.32583e-5, -2.50356, 2.79656, -1.79401, 0.714888, -0.165924, 0.0118379, 0.0037806, -1.57905e-4, 1.45323e-6, 0.11};

static const double pars_B97x[N_PARS] =
  {0.0027355, 0.0432970, -0.0669379, 0.0699060, -0.0474635, 0.0153092,
   15.8279, -26.8145, 17.8127, -5.98246, 1.25408, -0.270783, 0.0919536, -0.0140960, 0.0045466, 0.11};

static const double pars_cx13[N_PARS] =
  { 0.0024387, -0.0041526,  0.0025826,  0.0000012, -0.0007582,  0.0002764,
   -2.2030319,  2.1759315, -1.2997841,  0.5347267, -0.1588798,  0.0367329, -0.0077318,  0.0012667,  0.0000008, 0.106};

static const double pars_pw86[N_PARS] =
  { 0.0000006,  0.0402647, -0.0353219,  0.0116112, -0.0001555,  0.0000504,
   -1.8779594,  1.5198811, -0.5383109,  0.1352399, -0.0428465,  0.0117903,  0.0033791, -0.0000493,  0.0000071, 0.106};

static const double pars_b86r[N_PARS] =
  { 0.0045620, -0.0087000,  0.0073696, -0.0030244,  0.0003868,  0.0000944,
   -2.2089330,  2.1968353, -1.2662249,  0.4689964, -0.1165714,  0.0207188, -0.0029772,  0.0005982,  0.0000047, 0.106};

static void
gga_x_hjs_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hjs_params));

  xc_hyb_init_hybrid(p, 0.0);
  p->hyb_type[0] = XC_HYB_NONE;
}

#include "maple2c/gga_exc/gga_x_hjs.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_pbe = {
  XC_GGA_X_HJS_PBE,
  XC_EXCHANGE,
  "HJS screened exchange PBE version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-12,
  {N_PARS, names, desc, pars_PBE, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_pbe_sol = {
  XC_GGA_X_HJS_PBE_SOL,
  XC_EXCHANGE,
  "HJS screened exchange PBE_SOL version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  5e-12,
  {N_PARS, names, desc, pars_PBE_sol, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_b88 = {
  XC_GGA_X_HJS_B88,
  XC_EXCHANGE,
  "HJS screened exchange B88 version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-7, /* densities smaller than 1e-7 yield NaNs */
  {N_PARS, names, desc, pars_B88, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_b97x = {
  XC_GGA_X_HJS_B97X,
  XC_EXCHANGE,
  "HJS screened exchange B97x version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {N_PARS, names, desc, pars_B97x, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_cx13 = {
  XC_GGA_X_HJS_CX13,
  XC_EXCHANGE,
  "HJS screened exchange CX13 version",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_025902, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-11,
  {N_PARS, names, desc, pars_cx13, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_pw86 = {
  XC_GGA_X_HJS_PW86,
  XC_EXCHANGE,
  "HJS screened exchange PW86 version",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_025902, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-11,
  {N_PARS, names, desc, pars_b86r, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_b86r = {
  XC_GGA_X_HJS_B86R,
  XC_EXCHANGE,
  "HJS screened exchange B86r version",
  XC_FAMILY_GGA,
  {&xc_ref_Shukla2022_041003, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-11,
  {N_PARS, names, desc, pars_b86r, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, &work_gga, NULL
};
