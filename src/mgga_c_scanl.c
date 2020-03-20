/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_SCANL          702 /* SCAN correlation */
#define XC_MGGA_C_SCANL_RVV10    703 /* SCAN correlation + rVV10 correlation */
#define XC_MGGA_C_SCANL_VV10     704 /* SCAN correlation +  VV10 correlation */

#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_c_scanl.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl = {
  XC_MGGA_C_SCANL,
  XC_CORRELATION,
  "Deorbitalized SCAN (SCAN-L) correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161,&xc_ref_Sun2015_036402, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  NULL, NULL, 
  NULL, NULL, work_mgga,
};


static void
mgga_c_scan_rvv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_SCANL};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 15.7;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl_rvv10 = {
  XC_MGGA_C_SCANL_RVV10,
  XC_CORRELATION,
  "SCAN-L + rVV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Peng2016_041005, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  mgga_c_scan_rvv10_init, NULL,
  NULL, NULL, NULL
};

static void
mgga_c_scan_vv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_SCANL};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 14.0;
  p->nlc_C = 0.0093;
}

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_scanl_vv10 = {
  XC_MGGA_C_SCANL_VV10,
  XC_CORRELATION,
  "SCAN-L + VV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Mejia2017_052512, &xc_ref_Mejia2018_115161, &xc_ref_Brandenburg2016_115144, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_VV10 | MAPLE2C_FLAGS,
  1e-20,
  {0, NULL, NULL, NULL, NULL},
  mgga_c_scan_vv10_init, NULL,
  NULL, NULL, NULL
};
