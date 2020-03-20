/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_SCANL         700 /* Deorbitalized SCAN exchange */
#define XC_MGGA_X_REVSCANL      701 /* Deorbitalized revSCAN exchange */

typedef struct{
  double c1, c2, d, k1;
} mgga_x_scan_params;

static const mgga_x_scan_params par_scanl = {0.667, 0.8, 1.24, 0.065};
static const mgga_x_scan_params par_revscanl = {0.607, 0.7, 1.37, 0.065};

static void 
mgga_x_scanl_init(xc_func_type *p)
{
  mgga_x_scan_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_scan_params));
  params = (mgga_x_scan_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_SCANL:
    memcpy(params, &par_scanl, sizeof(mgga_x_scan_params));
    break;
  case XC_MGGA_X_REVSCANL:
    memcpy(params, &par_revscanl, sizeof(mgga_x_scan_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_scan\n");
    exit(1);
  }  
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_scanl.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_scanl = {
  XC_MGGA_X_SCANL,
  XC_EXCHANGE,
  "Deorbitalized SCAN (SCAN-L) exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Sun2015_036402, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_scanl_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_revscanl = {
  XC_MGGA_X_REVSCANL,
  XC_EXCHANGE,
  "Deorbitalized revised SCAN (revSCAN-L) exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Mezei2018_2469, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_scanl_init, NULL,
  NULL, NULL, work_mgga,
};

