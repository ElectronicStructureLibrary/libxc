/*
 Copyright (C) 2006-2016 M.A.L. Marques
                    2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_LDA_C_CHACHIYO_MOD  307   /* Chachiyo simple 2 parameter correlation with modified scaling   */
#define XC_LDA_C_KARASIEV_MOD  308   /* Karasiev reparameterization of Chachiyo with modified scaling   */

typedef struct {
  double ap, bp, af, bf;
} lda_c_chachiyo_mod_params;

static lda_c_chachiyo_mod_params par_chachiyo = {-0.01554535, 20.4562557, -0.007772675, 27.4203609};
static lda_c_chachiyo_mod_params par_karasiev = {-0.01554535, 21.7392245, -0.007772675, 28.3559732};

static void 
lda_c_chachiyo_mod_init(xc_func_type *p)
{
  lda_c_chachiyo_mod_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(lda_c_chachiyo_mod_params));
  params = (lda_c_chachiyo_mod_params *) (p->params);

  switch(p->info->number){
  case XC_LDA_C_CHACHIYO_MOD:
    memcpy(params, &par_chachiyo, sizeof(lda_c_chachiyo_mod_params));
    break;
  case XC_LDA_C_KARASIEV_MOD:
    memcpy(params, &par_karasiev, sizeof(lda_c_chachiyo_mod_params));
    break;
  default:
    fprintf(stderr, "Internal error in lda_c_chachiyo_mod\n");
    exit(1);
  }
}

#include "decl_lda.h"
#include "maple2c/lda_exc/lda_c_chachiyo_mod.c"
#include "work_lda.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_chachiyo_mod = {
  XC_LDA_C_CHACHIYO_MOD,
  XC_CORRELATION,
  "Chachiyo simple 2 parameter correlation with modified spin scaling",
  XC_FAMILY_LDA,
  {&xc_ref_Chachiyo2018_00712, &xc_ref_Chachiyo2016_021101, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  lda_c_chachiyo_mod_init, NULL,
  work_lda, NULL, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_lda_c_karasiev_mod = {
  XC_LDA_C_KARASIEV_MOD,
  XC_CORRELATION,
  "Karasiev reparameterization of Chachiyo",
  XC_FAMILY_LDA,
  {&xc_ref_Chachiyo2018_00712, &xc_ref_Karasiev2016_157101, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-24,
  {0, NULL, NULL, NULL, NULL},
  lda_c_chachiyo_mod_init, NULL,
  work_lda, NULL, NULL
};
