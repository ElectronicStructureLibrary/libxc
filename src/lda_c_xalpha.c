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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"


/************************************************************************
 Slater's Xalpha functional

    Exc = alpha Ex
************************************************************************/

/* This correlation functional, added to the exchange functional, produces
a total exchange and correlation functional, Exc, equal to 3/2 * alpha * Ex 
Setting alpha equal to one gives the *usual* Slater Xalpha functional,
whereas alpha equal to 2/3 just leaves the exchange functional unchanged */

#define XC_LDA_C_XALPHA  6   /* Slater's Xalpha              */

static void lda_c_xalpha(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc, FLOAT *kc)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;
  FLOAT a = 1.5*p->alpha - 1.0;
  int i;

  XC(lda)(p->lda_aux, rho, ec, vc, fc, NULL);

  if(ec != NULL)
    (*ec) *= a;

  if(vc != NULL)
    for(i=0; i<p->nspin; i++) vc[i] *= a;

  if(fc != NULL){
    int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++) fc[i] *= a;
  }
}

/* These prototypes are needed for the declaration of func_info_lda_c_xalpha */
void XC(lda_c_xalpha_init_default)(void *p_);
void XC(lda_c_xalpha_end)(void *p_);

const XC(func_info_type) XC(func_info_lda_c_xalpha) = {
  XC_LDA_C_XALPHA,
  XC_CORRELATION,
  "Slater's Xalpha",
  XC_FAMILY_LDA,
  NULL,
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  XC(lda_c_xalpha_init_default),  /* init */
  XC(lda_c_xalpha_end),           /* end  */
  lda_c_xalpha                   /* lda */
};

void XC(lda_c_xalpha_init)(XC(lda_type) *p, int nspin, int dim, FLOAT alpha)
{
  p->info = &XC(func_info_lda_c_xalpha);
  p->nspin = nspin;
  p->dim   = dim;
  p->alpha = alpha;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_x_init)(p->lda_aux, nspin, dim, XC_NON_RELATIVISTIC);
}

void XC(lda_c_xalpha_init_default)(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  XC(lda_c_xalpha_init)(p, p->nspin, 3, 1.0);
}

void XC(lda_c_xalpha_end)(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  free(p->lda_aux);
}
