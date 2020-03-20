/*
 Copyright (C) 2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_GGA_C_CHACHIYO  309   /* Chachiyo simple GGA correlation */

typedef struct {
  double ap, bp, af, bf, h;
} gga_c_chachiyo_params;

static gga_c_chachiyo_params par_chachiyo = {-0.01554535, 20.4562557, -0.007772675, 27.4203609, 0.06672632};

static void 
gga_c_chachiyo_init(xc_func_type *p)
{
  gga_c_chachiyo_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_c_chachiyo_params));
  params = (gga_c_chachiyo_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_C_CHACHIYO:
    memcpy(params, &par_chachiyo, sizeof(gga_c_chachiyo_params));
    break;
  default:
    fprintf(stderr, "Internal error in gga_c_chachiyo\n");
    exit(1);
  }
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_c_chachiyo.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_c_chachiyo = {
  XC_GGA_C_CHACHIYO,
  XC_CORRELATION,
  "Chachiyo simple GGA correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Chachiyo2018_00712, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  gga_c_chachiyo_init, NULL,
  NULL, work_gga, NULL
};
