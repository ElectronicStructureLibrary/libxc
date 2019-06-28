/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_RLDA          688 /* Reparametrized local-density approximation */

#include "maple2c/mgga_exc/mgga_x_rlda.c"
#include "work_mgga_new.c"

const xc_func_info_type xc_func_info_mgga_x_rlda = {
  XC_MGGA_X_RLDA,
  XC_EXCHANGE,
  "Reparametrized local-density approximation",
  XC_FAMILY_MGGA,
  {&xc_ref_Campi1978_263, &xc_ref_Koehl1996_835, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_I_HAVE_ALL,
  1.0e-23,
  0, NULL, NULL,
  NULL, NULL,
  NULL, NULL, work_mgga,
};
