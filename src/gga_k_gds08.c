/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_K_GDS08     591 /* Combined analytical theory with Monte Carlo sampling */
#define XC_GGA_K_GHDS10    592 /* As GDS08 but for an electron gas with spin */

static void
gga_k_gds08_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_VW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  xc_mix_init(p, 2, funcs_id, funcs_coef);  
}

const xc_func_info_type xc_func_info_gga_k_gds08 = {
  XC_GGA_K_GDS08,
  XC_KINETIC,
  "Combined analytical theory with Monte Carlo sampling",
  XC_FAMILY_GGA,
  {&xc_ref_Ghiringhelli2008_073104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_k_gds08_init, NULL,
  NULL, NULL, NULL
};

static void
gga_k_ghds10_init(xc_func_type *p)
{
  static int    funcs_id  [2] = {XC_GGA_K_TFVW, XC_LDA_K_GDS08_WORKER};
  static double funcs_coef[2] = {1.0, 1.0};

  static double par_k_gds08[] = {1.02, 0.163};
  
  xc_mix_init(p, 2, funcs_id, funcs_coef);

  xc_func_set_ext_params(p->func_aux[1], par_k_gds08);
}

const xc_func_info_type xc_func_info_gga_k_ghds10 = {
  XC_GGA_K_GHDS10,
  XC_KINETIC,
  "As GDS08 but for an electron gas with spin",
  XC_FAMILY_GGA,
  {&xc_ref_Ghiringhelli2010_014106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  gga_k_ghds10_init, NULL,
  NULL, NULL, NULL
};
