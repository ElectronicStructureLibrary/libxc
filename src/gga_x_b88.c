/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_B88          106 /* Becke 88 */
#define XC_GGA_X_OPTB88_VDW   139 /* Becke 88 reoptimized to be used with vdW functional of Dion et al */
#define XC_GGA_X_MB88         149 /* Modified Becke 88 for proton transfer */
#define XC_GGA_X_EB88         271 /* Non-empirical (excogitated) B88 functional of Becke and Elliott */
#define XC_GGA_X_B88M         570 /* Becke 88 reoptimized to be used with mgga_c_tau1 */

typedef struct{
  double beta, gamma;
} gga_x_b88_params;


static void 
gga_x_b88_init(xc_func_type *p)
{
  gga_x_b88_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_b88_params));
  params = (gga_x_b88_params *) (p->params);
  
  /* value of beta in standard Becke 88 functional */
  switch(p->info->number){
  case XC_GGA_X_B88:
    /* deffault values set by set_ext_params */
    break;
  case XC_GGA_X_OPTB88_VDW:
    params->beta  = 0.00336865923905927;
    params->gamma = 6.98131700797731;
    break;
  case XC_GGA_X_MB88:
    params->beta  = 0.0011;
    params->gamma = 6.0;
    break;
  case XC_GGA_X_EB88:
    params->beta  = 0.0050/M_CBRT2;
    params->gamma = 6.0;
    break;
  case XC_GGA_X_B88M:
    params->beta  = 0.0045;
    params->gamma = 6.0;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_b88\n");
    exit(1);
  }
}

static func_params_type ext_params[] = {
  {"_beta", 0.0042, "beta/X_FACTOR_C is the coefficient of the gradient expansion"},
  {"_gamma", 6.0, "gamma should be 6 to get the right asymptotics of Ex"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_b88_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b88_params *) (p->params);

  params->beta  = get_ext_param(p->info->ext_params, ext_params, 0);
  params->gamma = get_ext_param(p->info->ext_params, ext_params, 1);
}


#include "maple2c/gga_exc/gga_x_b88.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1988_3098, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  2, ext_params, set_ext_params,
  gga_x_b88_init, NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_optb88_vdw = {
  XC_GGA_X_OPTB88_VDW,
  XC_EXCHANGE,
  "opt-Becke 88 for vdW",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2010_022201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_mb88 = {
  XC_GGA_X_MB88,
  XC_EXCHANGE,
  "Modified Becke 88 for proton transfer",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2009_14415, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_eb88 = {
  XC_GGA_X_EB88,
  XC_EXCHANGE,
  "Non-empirical (excogitated) B88 functional of Becke and Elliott",
  XC_FAMILY_GGA,
  {&xc_ref_Elliott2009_1485, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init,  NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_b88m = {
  XC_GGA_X_B88M,
  XC_EXCHANGE,
  "Becke 88 reoptimized to be used with tau1",
  XC_FAMILY_GGA,
  {&xc_ref_Proynov2000_10013, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init,  NULL, 
  NULL, work_gga, NULL
};
