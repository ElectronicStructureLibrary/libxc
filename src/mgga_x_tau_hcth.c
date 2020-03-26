/*
 Copyright (C) 2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_TAU_HCTH        205 /* tau-HCTH from Boese and Handy */
#define XC_HYB_MGGA_X_BMK         279 /* Boese-Martin for kinetics     */
#define XC_HYB_MGGA_X_TAU_HCTH    282 /* Hybrid version of tau-HCTH    */

const double tHCTH_cx_local [4] = {1.10734, -1.0534, 6.3491, -2.5531};
const double tHCTH_cx_nlocal[4] = {0.00110, -0.3041, 6.9543, -0.7235};

const double BMK_cx_local [4] = { 0.474302, 2.77701, -11.4230, 11.7167};
const double BMK_cx_nlocal[4] = {-0.192212, 4.73936, -26.6188, 22.4891};

const double hyb_tHCTH_cx_local [4] = { 0.86735,  0.3008, 1.2208,   0.1574};
const double hyb_tHCTH_cx_nlocal[4] = {-0.00230, -0.2849, 5.4146, -10.909};

typedef struct{
  double cx_local[4];
  double cx_nlocal[4];
} mgga_x_tau_hcth_params;


static void 
mgga_x_tau_hcth_init(xc_func_type *p)
{
  mgga_x_tau_hcth_params *params;
  int ii;

  assert(p != NULL);
  assert(p->params == NULL);

  p->params = libxc_malloc(sizeof(mgga_x_tau_hcth_params));
  params = (mgga_x_tau_hcth_params *)(p->params);

  for(ii = 0; ii < 4; ii++){
    switch(p->info->number){
    case XC_MGGA_X_TAU_HCTH:
      params->cx_local[ii]  = tHCTH_cx_local[ii];
      params->cx_nlocal[ii] = tHCTH_cx_nlocal[ii];
      break;
    case XC_HYB_MGGA_X_BMK:
      p->cam_alpha = 0.42;
      params->cx_local[ii]  = BMK_cx_local[ii];
      params->cx_nlocal[ii] = BMK_cx_nlocal[ii];
    break;
    case XC_HYB_MGGA_X_TAU_HCTH:
      p->cam_alpha = 0.15;
      params->cx_local[ii]  = hyb_tHCTH_cx_local[ii];
      params->cx_nlocal[ii] = hyb_tHCTH_cx_nlocal[ii];
      break;
    default:
      fprintf(stderr, "Internal error in mgga_tau_hcth\n");
      exit(1);
      break;
    }
  }
  
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_tau_hcth.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_tau_hcth = {
  XC_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "tau-HCTH from Boese and Handy",
  XC_FAMILY_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_tau_hcth_init, NULL, 
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_bmk = {
  XC_HYB_MGGA_X_BMK,
  XC_EXCHANGE,
  "Boese-Martin for kinetics",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-15,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_tau_hcth_init, NULL, 
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_tau_hcth = {
  XC_HYB_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "Hybrid version of tau-HCTH",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_tau_hcth_init,  NULL, 
  NULL, NULL, work_mgga,
};
