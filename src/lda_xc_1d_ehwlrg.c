/*
 Copyright (C) 2006-2009 M.A.L. Marques

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

/* these functionals is for the soft-Coulomb interaction with alpha=1 */

#define XC_LDA_XC_1D_EHWLRG_1     536 /* LDA constructed from slab-like systems of 1 electron  */
#define XC_LDA_XC_1D_EHWLRG_2     537 /* LDA constructed from slab-like systems of 2 electrons */
#define XC_LDA_XC_1D_EHWLRG_3     538 /* LDA constructed from slab-like systems of 3 electrons */

typedef struct {
  FLOAT alpha;
  FLOAT a1, a2, a3;
} lda_xc_1d_ehwlrg_params;

static void 
lda_xc_1d_ehwlrg_init(XC(func_type) *p)
{
  lda_xc_1d_ehwlrg_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_xc_1d_ehwlrg_params));
  params = (lda_xc_1d_ehwlrg_params *) (p->params);

  switch(p->info->number){
  case XC_LDA_XC_1D_EHWLRG_1:
    params->alpha =  0.638;
    params->a1    = -0.803;
    params->a2    =  0.82;
    params->a3    = -0.47;
    break;
  case XC_LDA_XC_1D_EHWLRG_2:
    params->alpha =  0.604;
    params->a1    = -0.74;
    params->a2    =  0.68;
    params->a3    = -0.38;
    break;
  case XC_LDA_XC_1D_EHWLRG_3:
    params->alpha =  0.61;
    params->a1    = -0.77;
    params->a2    =  0.79;
    params->a3    = -0.48;
    break;
  }
}

#define XC_DIMENSIONS 1

#include "maple2c/lda_xc_1d_ehwlrg.c"

#define func maple2c_func
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_xc_1d_ehwlrg_1) = {
  XC_LDA_XC_1D_EHWLRG_1,
  XC_CORRELATION,
  "LDA constructed from slab-like systems of 1 electron",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D |  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_xc_1d_ehwlrg_2) = {
  XC_LDA_XC_1D_EHWLRG_2,
  XC_CORRELATION,
  "LDA constructed from slab-like systems of 2 electrons",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D |  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};


const XC(func_info_type) XC(func_info_lda_xc_1d_ehwlrg_3) = {
  XC_LDA_XC_1D_EHWLRG_3,
  XC_CORRELATION,
  "LDA constructed from slab-like systems of 3 electrons",
  XC_FAMILY_LDA,
  {&xc_ref_Entwistle2016_205134, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D |  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  lda_xc_1d_ehwlrg_init, NULL,
  work_lda, NULL, NULL
};
