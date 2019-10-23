/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_C_PBE_VWN      216 /* Perdew, Burke & Ernzerhof correlation based on VWN LDA */

typedef struct{
  double beta, gamma, BB;
} gga_c_pbe_vwn_params;

static void gga_c_pbe_vwn_init(xc_func_type *p)
{
  gga_c_pbe_vwn_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_c_pbe_vwn_params));
  params = (gga_c_pbe_vwn_params *) (p->params);

  /* most functionals have the same gamma and B */
  params->gamma = (1.0 - log(2.0))/(M_PI*M_PI);
  params->BB = 1.0; 

  switch(p->info->number){
  case XC_GGA_C_PBE_VWN:
    /* default set by set_ext_params */
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_pbe_vwn\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_beta",  0.06672455060314922,     "beta constant"},
  {"_gamma", 0.031090690869654895034, "(1 - ln(2))/Pi^2 in the PBE"},
  {"_B",     1.0, "Multiplies the A t^2 term. Used in the SPBE functional"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_c_pbe_vwn_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_c_pbe_vwn_params *) (p->params);

  params->beta  = get_ext_param(p->info->ext_params, ext_params, 0);
  params->gamma = get_ext_param(p->info->ext_params, ext_params, 1);
  params->BB    = get_ext_param(p->info->ext_params, ext_params, 2);
}

#include "maple2c/gga_exc/gga_c_pbe_vwn.c"
#include "work_gga.c"

const xc_func_info_type xc_func_info_gga_c_pbe_vwn = {
  XC_GGA_C_PBE_VWN,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof based on VWN correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Kraisler2010_042516, &xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1e-12,
  3, ext_params, set_ext_params,
  gga_c_pbe_vwn_init, NULL, 
  NULL, work_gga, NULL
};
