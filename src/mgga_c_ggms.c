/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_C_GGMS_PBE       736 /* localized-density mass approximation (PBE version) */
#define XC_MGGA_C_GGMS_PBE_SOL   737 /* localized-density mass approximation (PBEsol version) */

/* these parameters are simply copied from gga_c_pbe */
typedef struct{
  double beta, gamma, BB;
} mgga_c_ggms_params;

static void mgga_c_ggms_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_c_ggms_params));
}

#define PBE_N_PAR 3
static const char  *pbe_names[PBE_N_PAR]  = {"_beta", "_gamma", "_B"};
static const char  *pbe_desc[PBE_N_PAR]   = {
  "beta constant",
  "(1 - ln(2))/Pi^2 in the PBE",
  "Multiplies the A t^2 term. Used in the SPBE functional"};
static const double pbe_values[PBE_N_PAR] =
  {0.06672455060314922, 0.031090690869654895034, 1.0};
static const double pbe_sol_values[PBE_N_PAR] =
  {0.046, 0.031090690869654895034, 1.0};


#include "maple2c/mgga_exc/mgga_c_ggms.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_ggms_pbe = {
  XC_MGGA_C_GGMS_PBE,
  XC_CORRELATION,
  "localized-density mass approximation (PBE version)",
  XC_FAMILY_MGGA,
  {&xc_ref_Pitallis2022, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_values, set_ext_params_cpy},
  NULL, NULL,
  NULL, NULL, &work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_c_ggms_pbe_sol = {
  XC_MGGA_C_GGMS_PBE_SOL,
  XC_CORRELATION,
  "localized-density mass approximation (PBEsol version)",
  XC_FAMILY_MGGA,
  {&xc_ref_Pitallis2022, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {PBE_N_PAR, pbe_names, pbe_desc, pbe_sol_values, set_ext_params_cpy},
  NULL, NULL,
  NULL, NULL, &work_mgga,
};
