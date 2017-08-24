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

#define XC_LDA_C_HL   4   /* Hedin & Lundqvist            */
#define XC_LDA_C_GL   5   /* Gunnarson & Lundqvist        */
#define XC_LDA_C_vBH 17   /* von Barth & Hedin            */

typedef struct {
  FLOAT r[2], c[2];
} lda_c_hl_params;

static const lda_c_hl_params par_hl = { /* HL unpolarized only*/
  {21.0, 21.0}, {0.0225, 0.0225}
};

static const lda_c_hl_params par_gl = {
  {11.4, 15.9}, {0.0333, 0.0203}
};

static const lda_c_hl_params par_vbh = {
  {30.0, 75.0}, {0.0252, 0.0127}
};

static void 
lda_c_hl_init(XC(func_type) *p)
{
  lda_c_hl_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_c_hl_params));
  params = (lda_c_hl_params *) (p->params);

  switch(p->info->number){
  case XC_LDA_C_HL:
    memcpy(params, &par_hl, sizeof(lda_c_hl_params));
    break;
  case XC_LDA_C_GL:
    memcpy(params, &par_gl, sizeof(lda_c_hl_params));
    break;
  case XC_LDA_C_vBH:
    memcpy(params, &par_vbh, sizeof(lda_c_hl_params));
    break;
  default:
    fprintf(stderr, "Internal error in lda_c_hl\n");
    exit(1);
  }
}

#include "maple2c/lda_c_hl.c"

#define func maple2c_func
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_hl) = {
  XC_LDA_C_HL,
  XC_CORRELATION,
  "Hedin & Lundqvist",
  XC_FAMILY_LDA,
  {&xc_ref_Hedin1971_2064, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-16, 0.0,
  0, NULL, NULL,
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_c_gl) = {
  XC_LDA_C_GL,
  XC_CORRELATION,
  "Gunnarson & Lundqvist",
  XC_FAMILY_LDA,
  {&xc_ref_Gunnarsson1976_4274, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-12, 0.0,
  0, NULL, NULL,
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};

const XC(func_info_type) XC(func_info_lda_c_vbh) = {
  XC_LDA_C_vBH,
  XC_CORRELATION,
  "von Barth & Hedin",
  XC_FAMILY_LDA,
  {&xc_ref_vonBarth1972_1629, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-14, 0.0,
  0, NULL, NULL,
  lda_c_hl_init, NULL,
  work_lda, NULL, NULL
};
