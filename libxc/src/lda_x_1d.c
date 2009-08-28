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

#define XC_LDA_X_1D          21 /* Exchange in 1D     */

/* Warning: spin is missing! */

typedef struct{
  int interaction;  /* 0: exponentially screened; 1: soft-Coulomb */
  FLOAT bb;         /* screening parameter beta */
} lda_x_1d_params;

static void 
lda_x_1d_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(lda_x_1d_params));

  /* default value is soft-Coulomb with beta=1.0 */
  XC(lda_x_1d_set_params)(p, 1, 1.0);
}

static void 
lda_x_1d_end(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}

void 
XC(lda_x_1d_set_params)(XC(lda_type) *p, int interaction, FLOAT bb)
{
  lda_x_1d_params *params = (lda_x_1d_params *)(p->params);

  assert(params != NULL);

  params->interaction = interaction;
  params->bb          = bb;
}


static inline FLOAT FT_inter(FLOAT x, int interaction)
{
  if(interaction == 1)
    return 2.0*bessk0(x); 
  
  return 0.0;
}


static void func1(FLOAT *x, int n, void *ex)
{
  int interaction = *(int *)ex;
  int ii;
  
  for(ii=0; ii<n; ii++)
    x[ii] = FT_inter(x[ii], interaction);
}


static void func2(FLOAT *x, int n, void *ex)
{
  int interaction = *(int *)ex;
  int ii;
  
  for(ii=0; ii<n; ii++)
    x[ii] = x[ii]*FT_inter(x[ii], interaction);
}


static inline void
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  int interaction;
  FLOAT bb, bb2, R, int1, int2;

  assert(p->params != NULL);
  interaction = ((lda_x_1d_params *)p->params)->interaction;
  bb  =         ((lda_x_1d_params *)p->params)->bb;
  bb2 = bb*bb;

  R = M_PI*bb/(2.0*r->rs[1]);

  int1 = integrate(func1, (void *)(&interaction), 0.0, R);
  int2 = integrate(func2, (void *)(&interaction), 0.0, R);

  r->zk = -1.0/(2.0*M_PI*bb)*(int1 - int2/R);

  if(r->order < 1) return;

  r->dedrs = 1.0/(r->rs[2]*bb2) * FT_inter(R, interaction) * (1.0 - 1.0/bb2)
    + int2/(M_PI*M_PI*bb2);

  if(r->order < 2) return;

  /* TODO : second derivatives */
}

#define XC_DIMENSIONS 1
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_x_1d) = {
  XC_LDA_X_1D,
  XC_EXCHANGE,
  "Echange in 1D",
  XC_FAMILY_LDA,
  "Unpublished",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_x_1d_init,    /* init */
  lda_x_1d_end,     /* end  */
  work_lda,         /* lda  */
};
