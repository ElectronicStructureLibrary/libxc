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

typedef struct{
  double a[6], b[9]; /* pointers to the a and b parameters */
} gga_x_hjs_params;

static const double a_PBE[] = 
  {0.0159941, 0.0852995, -0.160368, 0.152645, -0.0971263, 0.0422061};
static const double b_PBE[] = 
  {5.33319, -12.4780, 11.0988, -5.11013, 1.71468, -0.610380, 0.307555, -0.0770547, 0.0334840};

static const double a_PBE_sol[] = 
  {0.0047333, 0.0403304, -0.0574615, 0.0435395, -0.0216251, 0.0063721};
static const double b_PBE_sol[] = 
  {8.52056, -13.9885, 9.28583, -3.27287, 0.843499, -0.235543, 0.0847074, -0.0171561, 0.0050552};

static const double a_B88[] =
  {0.00968615, -0.0242498, 0.0259009, -0.0136606, 0.00309606, -7.32583e-5};
static const double b_B88[] =
  {-2.50356, 2.79656, -1.79401, 0.714888, -0.165924, 0.0118379, 0.0037806, -1.57905e-4, 1.45323e-6};

static const double a_B97x[] =
  {0.0027355, 0.0432970, -0.0669379, 0.0699060, -0.0474635, 0.0153092};
static const double b_B97x[] =
  {15.8279, -26.8145, 17.8127, -5.98246, 1.25408, -0.270783, 0.0919536, -0.0140960, 0.0045466};

static const double a_B88_V2[] =
  {0.0253933, -0.0673075, 0.0891476, -0.0454168, -0.00765813, 0.0142506};
static const double b_B88_V2[] =
  {-2.6506, 3.91108, -3.31509, 1.54485, -0.198386, -0.136112, 0.0647862, 0.0159586, -0.000245066};

static void
gga_x_hjs_init(xc_func_type *p)
{
  gga_x_hjs_params *params;
  int ii;
  
  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_hjs_params));
  params = (gga_x_hjs_params *) (p->params);
  
  /* omega = 0.11 is set by ext_params */
  switch(p->info->number){
  case XC_GGA_X_HJS_PBE:
    for(ii = 0; ii < 6; ii++) params->a[ii] = a_PBE[ii];
    for(ii = 0; ii < 9; ii++) params->b[ii] = b_PBE[ii];
    break;
  case XC_GGA_X_HJS_PBE_SOL:
    for(ii = 0; ii < 6; ii++) params->a[ii] = a_PBE_sol[ii];
    for(ii = 0; ii < 9; ii++) params->b[ii] = b_PBE_sol[ii];
    break;
  case XC_GGA_X_HJS_B88:
    for(ii = 0; ii < 6; ii++) params->a[ii] = a_B88[ii];
    for(ii = 0; ii < 9; ii++) params->b[ii] = b_B88[ii];
    break;
  case XC_GGA_X_HJS_B97X:
    for(ii = 0; ii < 6; ii++) params->a[ii] = a_B97x[ii];
    for(ii = 0; ii < 9; ii++) params->b[ii] = b_B97x[ii];
    break;
  case XC_GGA_X_HJS_B88_V2:
    for(ii = 0; ii < 6; ii++) params->a[ii] = a_B88_V2[ii];
    for(ii = 0; ii < 9; ii++) params->b[ii] = b_B88_V2[ii];
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_hjs_init\n");
    exit(1);
  }
}

static const char  *omega_names[]  = {"omega"};
static const char  *omega_desc[]   = {"screening parameter"};
static const double omega_values[] = {0.11};

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
  {1, omega_names, omega_desc, omega_values, set_ext_params_omega},
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
  {1, omega_names, omega_desc, omega_values, set_ext_params_omega},
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
  {1, omega_names, omega_desc, omega_values, set_ext_params_omega},
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
  {1, omega_names, omega_desc, omega_values, set_ext_params_omega},
  gga_x_hjs_init, NULL, 
  NULL, work_gga, NULL
};
