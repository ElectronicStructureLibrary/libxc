/*
 Copyright (C) 2008 Georg Madsen
               2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_EV93     35 /* Engel and Vosko */
#define XC_GGA_X_ECMV92  215 /* Engel, Chevary, Macdonald, and Vosko */

typedef struct{
  /* numerator */
  double a1;
  double a2;
  double a3;
  /* denominator */
  double b1;
  double b2;
  double b3;
} gga_x_ev93_params;

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_ev93.c"
#include "work_gga.c"

static void
gga_x_ev93_init(xc_func_type *p)
{
  gga_x_ev93_params *params;

  assert(p != NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_ev93_params));
  params = (gga_x_ev93_params *) (p->params);

  switch(p->info->number) {
  case XC_GGA_X_EV93:
    /* default set by set_ext_params */
    break;
  case XC_GGA_X_ECMV92:
    params->a1=27.8428;
    params->a2=11.7683;
    params->a3=0.0;
    params->b1=27.5026;
    params->b2=5.7728;
    params->b3=0.0;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_ev93\n");
    exit(1);
  }
}

static const func_params_type ext_params[] = {
  {"_a1", 1.647127, "a1"},
  {"_a2", 0.980118, "a2"},
  {"_a3", 0.017399, "a3"},
  {"_b1", 1.523671, "a4"},
  {"_b2", 0.367229, "a5"},
  {"_b3", 0.011282, "a6"}
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_ev93_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ev93_params *) (p->params);

  params->a1 = get_ext_param(p->info->ext_params, ext_params, 0);
  params->a2 = get_ext_param(p->info->ext_params, ext_params, 1);
  params->a3 = get_ext_param(p->info->ext_params, ext_params, 2);
  params->b1 = get_ext_param(p->info->ext_params, ext_params, 3);
  params->b2 = get_ext_param(p->info->ext_params, ext_params, 4);
  params->b3 = get_ext_param(p->info->ext_params, ext_params, 5);
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ev93 = {
  XC_GGA_X_EV93,
  XC_EXCHANGE,
  "Engel and Vosko",
  XC_FAMILY_GGA,
  {&xc_ref_Engel1993_13164, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-25,
  6, ext_params, set_ext_params,
  gga_x_ev93_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_ecmv92 = {
  XC_GGA_X_ECMV92,
  XC_EXCHANGE,
  "Engel, Chevary, Macdonald and Vosko",
  XC_FAMILY_GGA,
  {&xc_ref_Engel1992_7, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-25,
  0, NULL, NULL,
  gga_x_ev93_init, NULL,
  NULL, work_gga, NULL
};
