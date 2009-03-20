/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_LDA_C_1D_CSC          18 /* Casula, Sorella, and Senatore 1D correlation     */

static void 
lda_c_1d_csc_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(int));

  /* default value is 1.0 for no particular reason */
  XC(lda_c_1d_csc_set_params)(p, 1.0);
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
XC(lda_c_1d_csc_set_params)(XC(lda_type) *p, FLOAT bb)
{
  int ii;

  assert(p->params != NULL);

  if     (bb == 0.1)
    ii = 0;
  else if(bb == 0.3)
    ii = 1;
  else if(bb == 0.5)
    ii = 2;
  else if(bb == 0.75)
    ii = 3;
  else if(bb == 1.0)
    ii = 4;
  else if(bb == 2.0)
    ii = 5;
  else if(bb == 4.0)
    ii = 6;
  else{
    fprintf(stderr, "Invalid value of parameter b = %f in lda_c_1d_csc_set_params", bb);
    exit(1);
  }

  *((int *)p->params) = ii;
}

static inline void
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  static const struct {
    FLOAT A, B, C, n, alpha, beta, m;
  } pp[] = {
    {  4.66,  2.092, 3.735, 1.379, 23.63,  109.9,    1.837},
    {  9.5,   1.85,  5.64,  0.882,  5.346,   6.69,   3.110},
    { 16.40,  2.90,  6.235, 0.908,  3.323,   2.23,   3.368},
    { 22.53,  2.09,  7.363, 0.906,  2.029,   0.394,  4.070},
    { 32.1,   3.77,  7.576, 0.941,  1.63,    0.198,  4.086},
    {110.5,   7.90,  8.37,  1.287,  1.399,   0.0481, 4.260},
    {413.0,  10.8,   7.99,  1.549,  1.308,   0.0120, 4.165}
  };

  int ii;
  FLOAT rs_n, rs_m, arg, larg, den;
  FLOAT darg, dden;

  assert(p->params != NULL);
  ii = *((int *)p->params);

  rs_n = POW(r->rs[1], pp[ii].n);
  rs_m = POW(r->rs[1], pp[ii].m);

  arg  = 1.0 + pp[ii].alpha*r->rs[1] + pp[ii].beta*rs_m;
  larg = LOG(arg);

  den  = pp[ii].A + pp[ii].B*rs_n + pp[ii].C*r->rs[2];

  r->zk  = -r->rs[1]*larg/den;
  r->zk /= 2.0; /* conversion from Ry to Hartree */

  if(r->order < 1) return;

  darg = pp[ii].alpha*r->rs[1] + pp[ii].beta*pp[ii].m*rs_m; /* times rs */
  dden = pp[ii].B*pp[ii].n*rs_n + 2.0*pp[ii].C*r->rs[2];    /* times rs */

  r->dedrs  = -((larg + darg/arg)*den - dden*larg)/(den*den);
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
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_c_1d_csc_init,    /* init */
  lda_c_1d_csc_end,     /* end  */
  work_lda,             /* lda  */
};
