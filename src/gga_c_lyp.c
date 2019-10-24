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
  gga_c_lyp_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_lyp_params));
  params = (gga_c_lyp_params *) (p->params);      

  /* values of constants in standard LYP functional */
  switch(p->info->number){
  case XC_GGA_C_LYP:
    /* default set by set_ext_params */
    break;
  case XC_GGA_C_TM_LYP:
    params->A = 0.0393;
    params->B = 0.21;
    params->c = 0.41;
    params->d = 0.15;
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_lyp\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_A", 0.04918, "Parameter A of LYP"},
  {"_B", 0.132,   "Parameter B of LYP"},
  {"_c", 0.2533,  "Parameter c of LYP"},
  {"_d", 0.349,   "Parameter d of LYP"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_c_lyp_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_c_lyp_params *) (p->params);

  params->A = get_ext_param(p->info->ext_params, ext_params, 0);
  params->B = get_ext_param(p->info->ext_params, ext_params, 1);
  params->c = get_ext_param(p->info->ext_params, ext_params, 2);
  params->d = get_ext_param(p->info->ext_params, ext_params, 3);
}

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
  4, ext_params, set_ext_params,
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
  0, NULL, NULL,
  xc_gga_c_lyp_init, NULL,
  NULL, work_gga, NULL
};
