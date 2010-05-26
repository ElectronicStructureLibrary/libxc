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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_LDA_C_1D_CSC          18 /* Casula, Sorella, and Senatore 1D correlation     */

typedef struct{
  int interaction;  /* 0: exponentially screened; 1: soft-Coulomb */
  int ii;           /* index in the parameter list */

  FLOAT bb;         /* screening parameter */
} lda_c_1d_csc_params;

static void 
lda_c_1d_csc_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(lda_c_1d_csc_params));

  /* default value is soft-Coulomb with beta=1.0 */
  XC(lda_c_1d_csc_set_params_)(p, 1, 1.0);
}

static void 
lda_c_1d_csc_end(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}

void 
XC(lda_c_1d_csc_set_params)(XC(func_type) *p, int interaction, FLOAT bb)
{
  assert(p != NULL && p->lda != NULL);
  XC(lda_c_1d_csc_set_params_)(p->lda, interaction, bb);
}

void 
XC(lda_c_1d_csc_set_params_)(XC(lda_type) *p, int interaction, FLOAT bb)
{
  lda_c_1d_csc_params *params = (lda_c_1d_csc_params *)(p->params);

  assert(params != NULL);

  params->interaction = -1;
  params->ii          = -1;

  if(interaction == 0){
    if     (bb == 0.1)
      params->ii = 0;
    else if(bb == 0.3)
      params->ii = 1;
    else if(bb == 0.5)
      params->ii = 2;
    else if(bb == 0.75)
      params->ii = 3;
    else if(bb == 1.0)
      params->ii = 4;
    else if(bb == 2.0)
      params->ii = 5;
    else if(bb == 4.0)
      params->ii = 6;
  }else if(interaction == 1){
    if     (bb == 0.5)
      params->ii = 7 + 0;
    if     (bb == 1.0)
      params->ii = 7 + 1;
  }

  if(params->ii < 0){
    fprintf(stderr, "Invalid value of parameters (inter,b) = (%d,%f) in lda_c_1d_csc_set_params", interaction, bb);
    exit(1);
  }

  params->interaction = interaction;
  params->bb          = bb;
}

static inline void
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  static const struct {
    FLOAT A, B, C, D, E, n1, n2, alpha, beta, m;
  } pp[] = {
    {  4.66,  0.0,  2.092, 3.735, 0.0, 1.379, 2.0, 23.63,  109.9,    1.837}, /* exponentially screened interaction */
    {  9.5,   0.0,  1.85,  5.64,  0.0, 0.882, 2.0,  5.346,   6.69,   3.110},
    { 16.40,  0.0,  2.90,  6.235, 0.0, 0.908, 2.0,  3.323,   2.23,   3.368},
    { 22.53,  0.0,  2.09,  7.363, 0.0, 0.906, 2.0,  2.029,   0.394,  4.070},
    { 32.1,   0.0,  3.77,  7.576, 0.0, 0.941, 2.0,  1.63,    0.198,  4.086},
    {110.5,   0.0,  7.90,  8.37,  0.0, 1.287, 2.0,  1.399,   0.0481, 4.260},
    {413.0,   0.0, 10.8,   7.99,  0.0, 1.549, 2.0,  1.308,   0.0120, 4.165},

    { 7.40, 1.120, 1.890, 0.0964,  0.0250,   2.0, 3.0, 2.431, 0.0142, 2.922}, /* soft-Coulomb interaction */
    {18.40, 0.0,   7.501, 0.10185, 0.012827, 2.0, 3.0, 1.511, 0.258,  4.424},
  };

  int ii;
  FLOAT rs_n1, rs_n2, rs_m, arg, larg, den, num;
  FLOAT darg, dden, dnum;

  assert(p->params != NULL);
  ii = ((lda_c_1d_csc_params *)p->params)->ii;

  rs_n1 = POW(r->rs[1], pp[ii].n1);
  rs_n2 = POW(r->rs[1], pp[ii].n2);
  rs_m  = POW(r->rs[1], pp[ii].m);

  arg  = 1.0 + pp[ii].alpha*r->rs[1] + pp[ii].beta*rs_m;
  larg = LOG(arg);

  den  = pp[ii].A + pp[ii].B*r->rs[1] + pp[ii].C*rs_n1 + pp[ii].D*rs_n2;
  num  = r->rs[1] + pp[ii].E*r->rs[2];

  r->zk  = -num*larg/den;
  r->zk /= 2.0; /* conversion from Ry to Hartree */

  if(r->order < 1) return;

  darg = pp[ii].alpha + pp[ii].beta*pp[ii].m*rs_m/r->rs[1];
  dden = pp[ii].B + pp[ii].C*pp[ii].n1*rs_n1/r->rs[1] + pp[ii].D*pp[ii].n2*rs_n2/r->rs[1];
  dnum = 2.0*pp[ii].E*r->rs[1];

  r->dedrs  = -((dnum*larg + num*darg/arg)*den - dden*num*larg)/(den*den);
  r->dedrs /= 2.0; /* conversion from Ry to Hartree */

  r->dedz   = 0.0; /* apparently the function is spin-unpolarized only */

  if(r->order < 2) return;

  /* TODO : second derivatives */
}

#define XC_DIMENSIONS 1
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_1d_csc) = {
  XC_LDA_C_1D_CSC,
  XC_CORRELATION,
  "Casula, Sorella & Senatore",
  XC_FAMILY_LDA,
  "M Casula, S Sorella, and G Senatore, Phys. Rev. B 74, 245427 (2006)",
  XC_FLAGS_1D |  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  lda_c_1d_csc_init,    /* init */
  lda_c_1d_csc_end,     /* end  */
  work_lda,             /* lda  */
};
