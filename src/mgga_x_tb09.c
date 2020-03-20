/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_BJ06         207 /* Becke & Johnson correction to Becke-Roussel 89  */
#define XC_MGGA_X_TB09         208 /* Tran & Blaha correction to Becke & Johnson  */
#define XC_MGGA_X_RPP09        209 /* Rasanen, Pittalis, and Proetto correction to Becke & Johnson  */

typedef struct{
  double c;
  double alpha;
} mgga_x_tb09_params;

static void
mgga_x_tb09_init(xc_func_type *p)
{
  mgga_x_tb09_params *params;

  p->params = libxc_malloc(sizeof(mgga_x_tb09_params));
  params = (mgga_x_tb09_params *)p->params;

  params->c = 0;

  switch(p->info->number){
  case XC_MGGA_X_BJ06:
    params->c = 1.0;
    params->alpha = 0.0;
    break;
  case XC_MGGA_X_TB09:
    /* the value of c should be passed by the calling code */
    params->alpha = 0.0;
    break;
  case XC_MGGA_X_RPP09:
    params->c = 1.0;
    params->alpha = 1.0;
    break;
  }
}

#define XC_NO_EXC
#include "decl_mgga.h"
#include "maple2c/mgga_vxc/mgga_x_tb09.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_bj06 = {
  XC_MGGA_X_BJ06,
  XC_EXCHANGE,
  "Becke & Johnson 06",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke2006_221101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};

static const func_params_type ext_params[] = {
  {"c", 1.0, "This parameter involves an average over the unit cell and must be calculated by the calling program."},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_tb09_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_tb09_params *) (p->params);

  params->c = get_ext_param(p->info->ext_params, ext_params, 0);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tb09 = {
  XC_MGGA_X_TB09,
  XC_EXCHANGE,
  "Tran & Blaha 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Tran2009_226401, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  1, ext_params, set_ext_params,
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_rpp09 = {
  XC_MGGA_X_RPP09,
  XC_EXCHANGE,
  "Rasanen, Pittalis & Proetto 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Rasanen2010_044112, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga,
};


