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

#define XC_LDA_X_1D          21 /* Exchange in 1D     */

typedef struct{
  int interaction;  /* 0: exponentially screened; 1: soft-Coulomb */
  FLOAT bb;         /* screening parameter beta */
} lda_x_1d_params;

static void 
lda_x_1d_init(XC(func_type) *p)
{
  assert(p->params == NULL);
  p->params = malloc(sizeof(lda_x_1d_params));
}


static inline FLOAT FT_inter(FLOAT x, int interaction)
{
  assert(interaction == 0 || interaction == 1);

  if(interaction == 0){
    FLOAT x2 = x*x;
    return expint_e1(x2)*EXP(x2);
  }else
    return 2.0*XC(bessel_K0)(x); 
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
func(const XC(func_type) *p, XC(lda_work_t) *r)
{
  static int spin_sign[2] = {+1, -1};
  static int spin_fact[2] = { 2,  1};

  int interaction, is;
  FLOAT bb, R, int1[2], int2[2];

  assert(p->params != NULL);
  interaction = ((lda_x_1d_params *)p->params)->interaction;
  bb  =         ((lda_x_1d_params *)p->params)->bb;

  r->f = 0.0;
  for(is=0; is<p->nspin; is++){
    R = M_PI*bb*(1.0 + spin_sign[is]*r->z)/(2.0*r->rs);

    if(R == 0.0) continue;

    int1[is] = XC(integrate)(func1, (void *)(&interaction), 0.0, R);
    int2[is] = XC(integrate)(func2, (void *)(&interaction), 0.0, R);

    r->f -= (1.0 + spin_sign[is]*r->z) *
      (int1[is] - int2[is]/R);
  }
  r->f *= spin_fact[p->nspin-1]/(4.0*M_PI*bb);

  if(r->order < 1) return;
  
  r->dfdrs = 0.0;
  r->dfdz  = 0.0;
  for(is=0; is<p->nspin; is++){
    if(1.0 + spin_sign[is]*r->z == 0.0) continue;

    r->dfdrs +=               int2[is];
    r->dfdz  -= spin_sign[is]*int1[is];
  }
  r->dfdrs *= spin_fact[p->nspin-1]/(2.0*M_PI*M_PI*bb*bb);
  r->dfdz  *= spin_fact[p->nspin-1]/(4.0*M_PI*bb);

  if(r->order < 2) return;

  r->d2fdrs2 = r->d2fdrsz = r->d2fdz2  = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT ft, aux = 1.0 + spin_sign[is]*r->z;

    if(aux == 0.0) continue;

    R  = M_PI*bb*aux/(2.0*r->rs);
    ft = FT_inter(R, interaction);
 
    r->d2fdrs2 -= aux*aux*ft;
    r->d2fdrsz += spin_sign[is]*aux*ft;
    r->d2fdz2  -= ft;
  }
  r->d2fdrs2 *= spin_fact[p->nspin-1]/(8.0*r->rs*r->rs*r->rs);
  r->d2fdrsz *= spin_fact[p->nspin-1]/(8.0*r->rs*r->rs);
  r->d2fdz2  *= spin_fact[p->nspin-1]/(8.0*r->rs);

  if(r->order < 3) return;

  /* TODO : third derivatives */
}

#define XC_DIMENSIONS 1
#include "work_lda.c"

static const func_params_type ext_params[] = {
  {  1, "Interaction: 0 (exponentially screened) | 1 (soft-Coulomb)"},
  {1.0, "Screening parameter beta"}
};

static void 
set_ext_params(XC(func_type) *p, const double *ext_params)
{
  lda_x_1d_params *params;
  double ff;

  assert(p != NULL && p->params != NULL);
  params = (lda_x_1d_params *)(p->params);

  ff = (ext_params == NULL) ? p->info->ext_params[0].value : ext_params[0];
  params->interaction = (int)round(ff);
  ff = (ext_params == NULL) ? p->info->ext_params[1].value : ext_params[1];
  params->bb = ff;

  assert(params->interaction == 0 || params->interaction == 1);
}

const XC(func_info_type) XC(func_info_lda_x_1d) = {
  XC_LDA_X_1D,
  XC_EXCHANGE,
  "Exchange in 1D",
  XC_FAMILY_LDA,
  {&xc_ref_Helbig2011_032503, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 0.0,
  2, ext_params, set_ext_params,
  lda_x_1d_init,    /* init */
  NULL,             /* end  */
  work_lda,         /* lda  */
  NULL,
  NULL
};
