/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_C_LYP    131  /* Lee, Yang & Parr */
#define XC_GGA_C_TM_LYP 559  /* Takkar and McCarthy reparametrization */

typedef struct{
  double A, B, c, d;
} gga_c_lyp_params;

void xc_gga_c_lyp_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_lyp_params));
}

#define LYP_N_PAR 4
static const char  *lyp_names[LYP_N_PAR]  = {"_A", "_B", "_c", "_d"};
static const char  *lyp_desc[LYP_N_PAR]   = {
  "Parameter A of LYP",
  "Parameter B of LYP",
  "Parameter c of LYP",
  "Parameter d of LYP"};
static const double lyp_values[LYP_N_PAR] =
  {0.04918, 0.132, 0.2533, 0.349};
static const double lyp_tm_values[LYP_N_PAR] =
  {0.0393, 0.21, 0.41, 0.15};
  
#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_lyp.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_lyp = {
  XC_GGA_C_LYP,
  XC_CORRELATION,
  "Lee, Yang & Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Lee1988_785, &xc_ref_Miehlich1989_200, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  {LYP_N_PAR, lyp_names, lyp_desc, lyp_values, set_ext_params_cpy},
  xc_gga_c_lyp_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_tm_lyp = {
  XC_GGA_C_TM_LYP,
  XC_CORRELATION,
  "Takkar and McCarthy reparametrization",
  XC_FAMILY_GGA,
  {&xc_ref_Thakkar2009_134109, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-32,
  {LYP_N_PAR, lyp_names, lyp_desc, lyp_tm_values, set_ext_params_cpy},
  xc_gga_c_lyp_init, NULL,
  NULL, work_gga, NULL
};
