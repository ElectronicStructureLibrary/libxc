/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_X_S12G                 495 /* Swart 2012 GGA exchange                 */
#define XC_HYB_GGA_X_S12H             496 /* Swart 2012 GGA hybrid exchange          */
#define XC_HYB_GGA_X_CAM_S12G         614 /* Swart 2012 range-separated GGA exchange */
#define XC_HYB_GGA_X_CAM_S12H         615 /* Swart 2012 range-separated GGA exchange */
	
typedef struct {
  double A, B, C, D, E;
  double bx;
} gga_x_s12_params;

static void
gga_x_s12_init(xc_func_type *p)
{
  gga_x_s12_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_s12_params));
  params = (gga_x_s12_params *) (p->params);

  params->bx  = 1.0; /* we initialize it here */

  switch(p->info->number){
  case XC_HYB_GGA_X_S12H:
    xc_hyb_init_hybrid(p, 0.0);
    break;
  case XC_HYB_GGA_X_CAM_S12G:
  case XC_HYB_GGA_X_CAM_S12H:
    xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
    break;
  }
}

#define S12G_N_PAR 5
static const char  *s12g_names[S12G_N_PAR]  = {"_A", "_B", "_C", "_D", "_E"};
static const char  *s12g_desc[S12G_N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter"
};
static const double s12g_values[S12G_N_PAR] = {
  1.03842032, 1.757-1.03842032, 0.00403198, 0.00104596, 0.00594635
};

#define S12H_N_PAR 6
static const char  *s12h_names[S12H_N_PAR]  = {"_A", "_B", "_C", "_D", "_E", "_alpha"};
static const char  *s12h_desc[S12H_N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter",
  "Fraction of exact exchange"
};
static const double s12h_values[S12H_N_PAR] = {
  1.02543951, 1.757-1.02543951, 0.00761554, 0.00211063, 0.00604672, 0.25
};

#define CAM_S12_N_PAR 8
static const char  *cam_s12_names[CAM_S12_N_PAR]  = {"_A", "_B", "_C", "_D", "_E", "_alpha", "_beta", "_omega"};
static const char  *cam_s12_desc[CAM_S12_N_PAR]   = {
  "A parameter",
  "B parameter",
  "C parameter",
  "D parameter",
  "E parameter",
  "fraction of HF exchange",
  "fraction of SR exchange",
  "range-separation parameter"
};

/* Paper reports alpha and beta in the Yanai convention, which differs from libxc */
static const double par_cam_s12g[CAM_S12_N_PAR] = {
  1.03323556, 1.757-1.03323556, 0.00417251, 0.00115216, 0.00706184, 0.34485046, -0.34485046, 1.52420731
};
static const double par_cam_s12h[CAM_S12_N_PAR] = {
  1.02149642, 1.757-1.02149642, 0.00825905, 0.00235804, 0.00654977, 0.25+0.10897845, -0.10897845, 0.48516891
};


static void
s12h_set_ext_params(xc_func_type *p, const double *ext_params)
{
  gga_x_s12_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_s12_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_S12G:
    set_ext_params_cpy(p, ext_params);
    break;
  case XC_HYB_GGA_X_S12H:
    set_ext_params_cpy_exx(p, ext_params);
    params->bx   = 1.0 - p->hyb_coeff[0];
    break;
  case XC_HYB_GGA_X_CAM_S12G:
  case XC_HYB_GGA_X_CAM_S12H:
    set_ext_params_cpy_cam(p, ext_params);
    params->bx   = 1.0 - (p->hyb_coeff[0] + p->hyb_coeff[1]);
    break;
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_x_s12.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_x_s12g = {
  XC_GGA_X_S12G,
  XC_EXCHANGE,
  "Swart 2012 GGA exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {S12G_N_PAR, s12g_names, s12g_desc, s12g_values, set_ext_params_cpy},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_s12h = {
  XC_HYB_GGA_X_S12H,
  XC_EXCHANGE,
  "Swart 2012 hybrid exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {S12H_N_PAR, s12h_names, s12h_desc, s12h_values, s12h_set_ext_params},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_cam_s12g = {
  XC_HYB_GGA_X_CAM_S12G,
  XC_EXCHANGE,
  "Swart 2012 range-separated hybrid GGA exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CAM_S12_N_PAR, cam_s12_names, cam_s12_desc, par_cam_s12g, s12h_set_ext_params},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_x_cam_s12h = {
  XC_HYB_GGA_X_CAM_S12H,
  XC_EXCHANGE,
  "Swart 2012 range-separated hybrid GGA exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2013_166, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-15,
  {CAM_S12_N_PAR, cam_s12_names, cam_s12_desc, par_cam_s12h, s12h_set_ext_params},
  gga_x_s12_init, NULL,
  NULL, work_gga, NULL
};
