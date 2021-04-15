/*
 Copyright (C) 2015 Narbe Mardirossian and Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_XC_WB97M_V   531 /* Mardirossian and Head-Gordon */

typedef struct{
  double omega;
} mgga_xc_wb97mv_params;


static void
hyb_mgga_xc_wb97mv_init(xc_func_type *p)
{
  mgga_xc_wb97mv_params *params;
  
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_xc_wb97mv_params));
  params = (mgga_xc_wb97mv_params *)(p->params);
  
  xc_hyb_init_cam(p, 1.0, -(1.0 - 0.15), 0.3);
  params->omega = p->hyb_params[1].sr.omega;
  
  p->hyb_number_terms = 3; /* we add a vv10 term */
  p->hyb_type[2] = XC_HYB_VDW_VV10;
  p->hyb_params[2].vv10.b = 6.0;
  p->hyb_params[2].vv10.C = 0.01;
}

#include "decl_mgga.h"
#include "maple2c/mgga_exc/hyb_mgga_xc_wb97mv.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_xc_wb97m_v = {
  XC_HYB_MGGA_XC_WB97M_V,
  XC_EXCHANGE_CORRELATION,
  "wB97M-V exchange-correlation functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Mardirossian2016_214110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-13,
  {0, NULL, NULL, NULL, NULL},
  hyb_mgga_xc_wb97mv_init, NULL,
  NULL, NULL, work_mgga,
};
