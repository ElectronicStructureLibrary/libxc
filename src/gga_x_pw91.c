/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_PW91         109 /* Perdew & Wang 91 */
#define XC_GGA_X_MPW91        119 /* Modified form of PW91 by Adamo & Barone */

typedef struct{
  double a, b, c, d, f, alpha, expo;
} gga_x_pw91_params;

static gga_x_pw91_params par_x_pw91 = /* b_PW91 ~ 0.0042 */
  {0.19645, 7.7956, 0.2743, -0.1508, 0.004, 100.0, 4.0};

static void 
gga_x_pw91_init(xc_func_type *p)
{
  gga_x_pw91_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pw91_params));
  params = (gga_x_pw91_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_PW91:
    memcpy(params, &par_x_pw91, sizeof(gga_x_pw91_params));
    break;
  case XC_GGA_X_MPW91:
    /* default set by set_ext_params */
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_pw91\n");
    exit(1);
  } 
}

/*
  === from nwchem source (xc_xmpw91.F) ===
  C. Adamo confirmed that there is a typo in the JCP paper
  b_mPW91 is 0.00426 instead of 0.0046
  also the power seems to be 3.72 and not 3.73
*/
static const func_params_type ext_params[] = {
  {"_bt",    0.00426, "a = 6 bt/X2S"},
  {"_alpha", 100.0,   "parameter of the exponential term"},
  {"_expo",  3.72,    "exponent of the power in the numerator"},
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_pw91_params *params;
  double bt, beta;
  
  assert(p != NULL && p->params != NULL);
  params = (gga_x_pw91_params *) (p->params);

  bt = get_ext_param(p->info->ext_params, ext_params, 0);
  params->alpha = get_ext_param(p->info->ext_params, ext_params, 1);
  params->expo  = get_ext_param(p->info->ext_params, ext_params, 2);


  beta         =  5.0*pow(36.0*M_PI,-5.0/3.0);
  params->a    =  6.0*bt/X2S;
  params->b    =  1.0/X2S;
  params->c    =  bt/(X_FACTOR_C*X2S*X2S);
  params->d    = -(bt - beta)/(X_FACTOR_C*X2S*X2S);
  params->f    =  1.0e-6/(X_FACTOR_C*pow(X2S, params->expo));
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_pw91.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_pw91 = {
  XC_GGA_X_PW91,
  XC_EXCHANGE,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1991, &xc_ref_Perdew1992_6671, &xc_ref_Perdew1992_6671_err, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  0, NULL, NULL,
  gga_x_pw91_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_mpw91 = {
  XC_GGA_X_MPW91,
  XC_EXCHANGE,
  "mPW91 of Adamo & Barone",
  XC_FAMILY_GGA,
  {&xc_ref_Adamo1998_664, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-31,
  3, ext_params, set_ext_params,
  gga_x_pw91_init, NULL,
  NULL, work_gga, NULL
};
