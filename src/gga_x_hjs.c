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
#define XC_GGA_X_HJS_APBE    648 /* HJS screened exchange APBE version */
#define XC_GGA_X_HJS_PBEINT  649 /* HJS screened exchange PBEint version */
#define XC_GGA_X_HJS_SG4     650 /* HJS screened exchange SG4 version */

typedef struct{
  double a[6], b[9]; /* pointers to the a and b parameters */
} gga_x_hjs_params;

#define N_PARS 16
static const char  *names[N_PARS]  = {"_a0", "_a1", "_a2", "_a3", "_a4", "_a5", "_b0", "_b1", "_b2", "_b3", "_b4", "_b5", "_b6", "_b7", "_b8", "_omega"};
static const char  *desc[N_PARS]   = {"a0", "a1", "a2", "a3", "a4", "a5", "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "omega"};

static const double pars_PBE[N_PARS] =
  {0.0159941, 0.0852995, -0.160368, 0.152645, -0.0971263, 0.0422061,
   5.33319, -12.4780, 11.0988, -5.11013, 1.71468, -0.610380, 0.307555, -0.0770547, 0.0334840,
   0.11};

static const double pars_PBE_sol[N_PARS] =
  {0.0047333, 0.0403304, -0.0574615, 0.0435395, -0.0216251, 0.0063721,
   8.52056, -13.9885, 9.28583, -3.27287, 0.843499, -0.235543, 0.0847074, -0.0171561, 0.0050552,
   0.11};

static const double pars_B88[N_PARS] =
  {0.00968615, -0.0242498, 0.0259009, -0.0136606, 0.00309606, -7.32583e-5,
   -2.50356, 2.79656, -1.79401, 0.714888, -0.165924, 0.0118379, 0.0037806, -1.57905e-4, 1.45323e-6,
   0.11};

static const double pars_B97x[N_PARS] =
  {0.0027355, 0.0432970, -0.0669379, 0.0699060, -0.0474635, 0.0153092,
   15.8279, -26.8145, 17.8127, -5.98246, 1.25408, -0.270783, 0.0919536, -0.0140960, 0.0045466,
   0.11};

static const double pars_APBE[N_PARS] =
  {0.02273221, -0.06039976, 0.07814309, -0.05481418, 0.01738706, 0.0008502,
   -2.91638499, 3.86022942, -2.97137272, 1.52515139, -0.57570311, 0.18417949, -0.0462847, 0.01406773, 0.00066305,
   0.37
  };

static const double pars_PBEINT[N_PARS] =
  {0.0047333, -0.0101441, 0.01072278, -0.00608029, 0.00144668, 0.0000766,
   -2.39737305, 2.5499878, -1.54388838, 0.59184499, -0.15548226, 0.03286721, -0.00608099, 0.00122022, 0.00005888,
    0.53
  };

static const double pars_SG4[N_PARS] =
  {0.02273221, 0.04855297, -0.12481836, 0.09963666, -0.03066345, 0.00432085,
   1.93055176, -5.47264766, 5.37181798, -2.91240963, 1.12772142, -0.31905254, 0.10263425, -0.02540923, 0.00343592,
    0.53
  };

static void
gga_x_hjs_init(xc_func_type *p)
{
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hjs_params));

  xc_hyb_init_hybrid(p, 0.0);
  p->hyb_type[0] = XC_HYB_NONE;
}

#include "decl_gga.h"
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
  NULL, work_gga, NULL
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
  NULL, work_gga, NULL
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
  NULL,  work_gga, NULL
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
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_apbe = {
  XC_GGA_X_HJS_APBE,
  XC_EXCHANGE,
  "HJS screened exchange APBE version",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {N_PARS, names, desc, pars_APBE, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_pbeint = {
  XC_GGA_X_HJS_PBEINT,
  XC_EXCHANGE,
  "HJS screened exchange PBEint version",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {N_PARS, names, desc, pars_PBEINT, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_hjs_sg4 = {
  XC_GGA_X_HJS_SG4,
  XC_EXCHANGE,
  "HJS screened exchange SG4 version",
  XC_FAMILY_GGA,
  {&xc_ref_Jana2019_042515, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-10,
  {N_PARS, names, desc, pars_SG4, set_ext_params_cpy_omega},
  gga_x_hjs_init, NULL,
  NULL, work_gga, NULL
};
