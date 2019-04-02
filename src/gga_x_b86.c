/*
 Copyright (C) 2006-2014 L. Talirz, M.A.L. Marques
                    2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_B86          103 /* Becke 86 Xalpha,beta,gamma                      */
#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_B86_R         41 /* Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_OPTB86B_VDW  171 /* Becke 86 reoptimized for use with vdW functional of Dion et al */

typedef struct{
  double beta, gamma, omega;
} gga_x_b86_params;


static void
gga_x_b86_init(xc_func_type *p)
{
  gga_x_b86_params *params;
  double mu, kappa;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_b86_params));
  params = (gga_x_b86_params *) (p->params);

  /* value of beta and gamma in Becke 86 functional */
  switch(p->info->number){
  case XC_GGA_X_B86:
    /* default set by set_ext_params */
    break;
  case XC_GGA_X_B86_MGC:
    params->beta  = 0.00375/X_FACTOR_C;
    params->gamma = 0.007;
    params->omega = 4.0/5.0;
    break;
  case XC_GGA_X_B86_R:
    mu = MU_GE;
    kappa = 0.7114;

    params->beta  = mu*X2S*X2S;
    params->gamma = mu*X2S*X2S/kappa;
    params->omega = 4.0/5.0;
    break;
  case XC_GGA_X_OPTB86B_VDW:
    mu = MU_GE;
    params->beta  = mu*X2S*X2S;
    params->gamma = mu*X2S*X2S;
    params->omega = 4.0/5.0;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_b86\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_beta", 0.0036/X_FACTOR_C, "Small x limit"},
  {"_gamma", 0.004, "Parameter in the denominator"},
  {"_omega", 1.0, "Exponent of denominator"},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_b86_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b86_params *) (p->params);

  params->beta  = get_ext_param(p->info->ext_params, ext_params, 0);
  params->gamma = get_ext_param(p->info->ext_params, ext_params, 1);
  params->omega = get_ext_param(p->info->ext_params, ext_params, 2);
}


#include "maple2c/gga_exc/gga_x_b86.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_b86 = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  3, ext_params, set_ext_params,
  gga_x_b86_init, NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_b86_mgc = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_x_b86_init, NULL, 
  NULL, work_gga, NULL
};

const xc_func_info_type xc_func_info_gga_x_b86_r = {
  XC_GGA_X_B86_R,
  XC_EXCHANGE,
  "Revised Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  {&xc_ref_Hamada2014_121103, &xc_ref_Becke1986_4524, &xc_ref_Becke1986_7184, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_x_b86_init, NULL, 
  NULL, work_gga, NULL

};

const xc_func_info_type xc_func_info_gga_x_optb86b_vdw = {
  XC_GGA_X_OPTB86B_VDW,
  XC_EXCHANGE,
  "Becke 86 reoptimized for use with vdW functional of Dion et al",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2011_195131, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_x_b86_init, NULL,
  NULL, work_gga, NULL
};
