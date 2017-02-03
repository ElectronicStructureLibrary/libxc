/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include "util.h"

/************************************************************************
 Correlation energy per particle and potential of a HEG as parametrized 
 by 
   Perdew & Zunger
   Ortiz & Ballone
************************************************************************/

#define XC_LDA_C_PZ       9   /* Perdew & Zunger              */
#define XC_LDA_C_PZ_MOD  10   /* Perdew & Zunger (Modified)   */
#define XC_LDA_C_OB_PZ   11   /* Ortiz & Ballone (PZ)         */

typedef struct {
  FLOAT gamma[2];
  FLOAT beta1[2];
  FLOAT beta2[2];
  FLOAT a[2], b[2], c[2], d[2];
} lda_c_pz_params;

static lda_c_pz_params pz_original = {
  {-0.1423, -0.0843},  /* gamma */
  { 1.0529,  1.3981},  /* beta1 */
  { 0.3334,  0.2611},  /* beta2 */
  { 0.0311,  0.01555}, /*  a    */
  {-0.048,  -0.0269},  /*  b    */
  { 0.0020,  0.0007},  /*  c    */
  {-0.0116, -0.0048}   /*  d    */
};

static lda_c_pz_params pz_modified = {
  {-0.1423, -0.0843},   
  { 1.0529,  1.3981}, 
  { 0.3334,  0.2611}, 
  { 0.0311,  0.01555},
  {-0.048,  -0.0269},   
  { 0.0020191519406228,  0.00069255121311694},
  {-0.0116320663789130, -0.00480126353790614}
};

static lda_c_pz_params pz_ob = {
  {-0.103756, -0.065951},
  { 0.56371,   1.11846},
  { 0.27358,   0.18797},
  { 0.031091,  0.015545},
  {-0.046644, -0.025599},
  { 0.00419,   0.00329},  /* the sign of c[0] and c[1] is different from [2], but is consistent
                             with the continuity requirement. There is nothing in [3] about this. */
  {-0.00983,  -0.00300}
};

static void 
lda_c_pz_init(XC(func_type) *p)
{
  lda_c_pz_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_c_pz_params));
  params = (lda_c_pz_params *) (p->params);

  switch(p->info->number){
  case XC_LDA_C_PZ:
    memcpy(params, &pz_original, sizeof(lda_c_pz_params));
    break;
  case XC_LDA_C_PZ_MOD:
    memcpy(params, &pz_modified, sizeof(lda_c_pz_params));
    break;
  case XC_LDA_C_OB_PZ:
    memcpy(params, &pz_ob, sizeof(lda_c_pz_params));
    break;
  default:
    fprintf(stderr, "Internal error in lda_c_pz\n");
    exit(1);
  }
}

#include "maple2c/lda_c_pz.c"

#define func maple2c_func
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_pz) = {
  XC_LDA_C_PZ,
  XC_CORRELATION,
  "Perdew & Zunger",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1981_5048, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_c_pz_mod) = {
  XC_LDA_C_PZ_MOD,
  XC_CORRELATION,
  "Perdew & Zunger (Modified)",
  XC_FAMILY_LDA,
  {&xc_ref_Perdew1981_5048_mod, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_c_ob_pz) = {
  XC_LDA_C_OB_PZ,
  XC_CORRELATION,
  "Ortiz & Ballone (PZ parametrization)",
  XC_FAMILY_LDA,
  {&xc_ref_Ortiz1994_1391, &xc_ref_Ortiz1994_1391_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_c_pz_init, NULL,
  work_lda, NULL, NULL
};
