/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_OL2          183 /* Exchange form based on Ou-Yang and Levy v.2 */

typedef struct{
  double aa, bb, cc;
} gga_x_ol2_params;

static void 
gga_x_ol2_init(xc_func_type *p)
{
  gga_x_ol2_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_ol2_params));
  params = (gga_x_ol2_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_OL2:
    params->aa = M_CBRT2*0.07064/X_FACTOR_C;
    params->bb = M_CBRT2*0.07064/X_FACTOR_C;
    params->cc = M_CBRT2*M_CBRT2*0.07064*34.0135/X_FACTOR_C;
    break;
  }
}

#include "maple2c/gga_exc/gga_x_ol2.c"
#include "work_gga_new.c"

const xc_func_info_type xc_func_info_gga_x_ol2 = {
  XC_GGA_X_OL2,
  XC_EXCHANGE,
  "Exchange form based on Ou-Yang and Levy v.2",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, &xc_ref_OuYang1991_379, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  5e-26,
  0, NULL, NULL,
  gga_x_ol2_init, NULL, 
  NULL, work_gga, NULL
};
