/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_RSCANL         705 /* Regularized SCAN exchange */

typedef struct{
  double c1, c2, d, k1;
  double taur, alphar; 
} mgga_x_rscan_params;

static const mgga_x_rscan_params par_rscan = {0.667, 0.8, 1.24, 0.065, 1.0e-4, 1.0e-3};

static void 
mgga_x_rscan_init(xc_func_type *p)
{
  mgga_x_rscan_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_x_rscan_params));
  params = (mgga_x_rscan_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_X_RSCANL:
    memcpy(params, &par_rscan, sizeof(mgga_x_rscan_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_rscan\n");
    exit(1);
  }  
}

#include "mgga_k_pcopt.c"
#include "maple2c/mgga_exc/mgga_x_rscan.c"
#include "work_mgga_deorb.c"

const xc_func_info_type xc_func_info_mgga_x_rscanl = {
  XC_MGGA_X_RSCANL,
  XC_EXCHANGE,
  "Deorbitalized and regularized SCAN (rSCAN-L) exchange",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Bartok2019_161101, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_HAVE_EXC,
  1e-12,
  0, NULL, NULL,
  mgga_x_rscan_init, NULL,
  NULL, NULL, work_mgga_deorb,
};
