/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_X_N12          82 /* N12 functional from Minnesota    */
#define XC_HYB_GGA_X_N12_SX   81 /* N12-SX functional from Minnesota */
#define XC_GGA_X_GAM          32 /* GAM functional from Minnesota */

static const double CC_N12[4][4] = {
  { 1.00000e+00,  5.07880e-01,  1.68233e-01,  1.28887e-01},
  { 8.60211e-02, -1.71008e+01,  6.50814e+01, -7.01726e+01},
  {-3.90755e-01,  5.13392e+01, -1.66220e+02,  1.42738e+02},
  { 4.03611e-01, -3.44631e+01,  7.61661e+01, -2.41834e+00}
};

static const double CC_N12_SX[4][4] = {
  /* Indices are wrong in the original paper; the first two indices
     need to be flipped */
  { 6.81116e-01,  1.88858e+00,  1.78590e+00,  8.79456e-01},
  {-8.12270e-02, -1.08723e+00, -4.18682e+00, -3.00000e+01},
  { 5.36236e-01, -5.45678e+00,  3.00000e+01,  5.51105e+01},
  {-7.09913e-01,  1.30001e+01, -7.24877e+01,  2.98363e+01}
};

static const double CC_GAM[4][4] = {
  { 1.32730,    0.886102, -5.73833,   8.60197},
  {-0.786018,  -4.78787,   3.90989,  -2.11611},
  { 0.802575,  14.4363,    8.42735,  -6.21552},
  {-0.142331, -13.4598,    1.52355, -10.0530}
};

typedef struct{
  double CC[4][4];
} gga_x_n12_params;


static void
gga_x_n12_init(xc_func_type *p)
{
  gga_x_n12_params *params;
  int ii, jj;

  assert(p != NULL);

  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_n12_params));
  params = (gga_x_n12_params *) (p->params);

  const double (*pCC)[4];
  
  switch(p->info->number){
  case XC_GGA_X_N12: 
    pCC = CC_N12;
    break;
  case XC_HYB_GGA_X_N12_SX:
    pCC = CC_N12_SX;
    p->cam_alpha = 0.00;
    p->cam_beta  = 0.25;
    p->cam_omega = 0.11;
    break;
  case XC_GGA_X_GAM:
    pCC = CC_GAM;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_n12\n");
    exit(1);
  }

  for(ii = 0; ii < 4; ii++){
    for(jj = 0; jj < 4; jj++){
      params->CC[ii][jj] = pCC[ii][jj];
    }
  }
  
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_n12.c"
#include "work_gga.c"


#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_n12 = {
  XC_GGA_X_N12,
  XC_EXCHANGE,
  "Minnesota N12 exchange functional",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2012_2310, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_n12_sx = {
  XC_HYB_GGA_X_N12_SX,
  XC_EXCHANGE,
  "Minnesota N12-SX exchange functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Peverati2012_16187, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_gam = {
  XC_GGA_X_GAM,
  XC_EXCHANGE,
  "Minnesota GAM exhange functional",
  XC_FAMILY_GGA,
  {&xc_ref_Yu2015_12146, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  gga_x_n12_init, NULL,
  NULL, work_gga, NULL
};
