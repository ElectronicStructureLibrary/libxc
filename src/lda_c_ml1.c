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

#define XC_LDA_C_ML1    22   /* Modified LSD (version 1) of Proynov and Salahub */
#define XC_LDA_C_ML2    23   /* Modified LSD (version 2) of Proynov and Salahub */

typedef struct {
  FLOAT fc, q;
} lda_c_ml1_params;

static void 
lda_c_ml1_init(XC(func_type) *p)
{
  lda_c_ml1_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_c_ml1_params));
  params = (lda_c_ml1_params *) (p->params);

  switch(p->info->number){
  case XC_LDA_C_ML1:
    params->fc = 0.2026;
    params->q  = 0.084;
    break;
  case XC_LDA_C_ML2:
    params->fc = 0.266;
    params->q  = 0.5;
    break;
  default:
    fprintf(stderr, "Internal error in lda_c_ml1\n");
    exit(1);
  }
}

#include "maple2c/lda_c_ml1.c"

#define func maple2c_func
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_ml1) = {
  XC_LDA_C_ML1,
  XC_CORRELATION,
  "Modified LSD (version 1) of Proynov and Salahub",
  XC_FAMILY_LDA,
  {&xc_ref_Proynov1994_7874, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  lda_c_ml1_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_c_ml2) = {
  XC_LDA_C_ML2,
  XC_CORRELATION,
  "Modified LSD (version 2) of Proynov and Salahub",
  XC_FAMILY_LDA,
  {&xc_ref_Proynov1994_7874, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  lda_c_ml1_init, NULL,
  work_lda, NULL, NULL
};
