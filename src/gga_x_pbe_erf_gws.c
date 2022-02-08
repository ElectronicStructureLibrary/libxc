/*
 2018 Authored by Andrea Kreppel
 2022 Edited by Henryk Laqua

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.

 Short-range PBE exchange functional Goll/Werner/Stoll
 E. Goll, H.-J. Werner, and H. Stoll., Phys. Chem. Chem. Phys. 7, 3917 (2005).
 DOI:10.1039/B509242F
*/

#include "util.h"

#define  XC_GGA_X_PBE_ERF_GWS                   742 /* Short ranged PBE exchange (erfc) */

#define N_PAR 4

typedef struct{
  double kappa, b_PBE, ax, omega;
} gga_x_pbe_erf_gws_params;


static const char  *param_names[N_PAR]  = {"_kappa", "_b_PBE", "_ax","_omega"};
static const char  *param_desc[N_PAR]   = {
  "kappa from PBE",
  "original b (AKA mu) from PBE",
  "exponent in short-range b-expansion",
  "range-separation screening parameter (AKA mu)",
};
static const double param_values[N_PAR] =
  {0.8040, 0.2195149727645171, 19.0, 0.5};


static void xc_gga_x_pbe_erf_gws_init(xc_func_type *p)
{
  xc_hyb_init_sr(p, 0.0,0.0);
  p->hyb_type[0] = XC_HYB_NONE;
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_erf_gws_params));
}

#include "maple2c/gga_exc/gga_x_pbe_erf_gws.c"
#include "work_gga.c"

const xc_func_info_type xc_func_info_gga_x_pbe_erf_gws = {
  XC_GGA_X_PBE_ERF_GWS,
  XC_EXCHANGE,
  "Short ranged PBE exchange (erfc)",
  XC_FAMILY_GGA,
  {&xc_ref_Goll2005_3917,&xc_ref_Goll2006_276,NULL,NULL,NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-13,
  {N_PAR, param_names, param_desc, param_values, set_ext_params_cpy_omega},
  xc_gga_x_pbe_erf_gws_init, NULL,
  NULL, work_gga, NULL
};

static void xc_hyb_gga_x_pbe_erf_gws_init(xc_func_type *p)
{
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_x_pbe_erf_gws_params));
}

const xc_func_info_type xc_func_info_hyb_gga_x_pbe_erf_gws = {
  XC_HYB_GGA_X_PBE_ERF_GWS,
  XC_EXCHANGE,
  "Hybrid Exchange functional with short-range PBE (GWS) exchange (erf)",
  XC_FAMILY_GGA,
  {&xc_ref_Goll2005_3917,&xc_ref_Goll2006_276,NULL,NULL,NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-13,
  {N_PAR, param_names, param_desc, param_values, set_ext_params_cpy_lc},
  xc_hyb_gga_x_pbe_erf_gws_init, NULL,
  NULL, work_gga, NULL
};

