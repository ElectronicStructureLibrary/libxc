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

#define XC_LDA_X         1   /* Exchange                     */
#define XC_LDA_C_XALPHA  6   /* Slater Xalpha                */

/*  
    Slater's Xalpha functional (Exc = alpha Ex)
    
    Note: this is to be added to the exchange

    This correlation functional, added to the exchange functional, produces
    a total exchange and correlation functional, Exc, equal to 3/2 * alpha * Ex 
    Setting alpha equal to one gives the *usual* Slater Xalpha functional,
    whereas alpha equal to 2/3 just leaves the exchange functional unchanged 
*/

static void 
lda_x_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(float));

  /* exchange is equal to xalpha with a parameter of 4/3 */
  XC(lda_c_xalpha_set_params)(p, 4.0/3.0);
}

static void 
lda_c_xalpha_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(float));

  /* This gives the usual Xalpha functional */
  XC(lda_c_xalpha_set_params)(p, 1.0);
}

static void 
lda_x_end(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}


void 
XC(lda_c_xalpha_set_params)(XC(lda_type) *p, FLOAT alpha)
{
  assert(p->params != NULL);

  *((FLOAT *)p->params) = 1.5*alpha - 1.0;
}


static inline void 
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  FLOAT alpha, ax, fz, dfz, d2fz, d3fz;

  assert(p->params != NULL);
  alpha = *((FLOAT *)p->params);
  
  ax = -alpha*0.458165293283142893475554485052; /* -alpha * 3/4*POW(3/(2*M_PI), 2/3) */
  
  r->zk = ax/r->rs[1];
  if(p->nspin == XC_POLARIZED){
    fz  = 0.5*(pow(1.0 + r->zeta,  4.0/3.0) + pow(1.0 - r->zeta,  4.0/3.0));
    r->zk *= fz;
  }

  if(r->order < 1) return;
  
  r->dedrs = -ax/r->rs[2];
  if(p->nspin == XC_POLARIZED){
    dfz = 2.0/3.0*(pow(1.0 + r->zeta,  1.0/3.0) - pow(1.0 - r->zeta,  1.0/3.0));

    r->dedrs *= fz;
    r->dedz   = ax/r->rs[1]*dfz;
  }

  if(r->order < 2) return;
    
  r->d2edrs2 = 2.0*ax/(r->rs[1]*r->rs[2]);
  if(p->nspin == XC_POLARIZED){
    if(ABS(r->zeta) == 1.0)
      d2fz = FLT_MAX;
    else
      d2fz = 2.0/9.0*(pow(1.0 + r->zeta,  -2.0/3.0) + pow(1.0 - r->zeta,  -2.0/3.0));
    
    r->d2edrs2 *= fz;
    r->d2edrsz = -ax/r->rs[2]*dfz;
    r->d2edz2  =  ax/r->rs[1]*d2fz;
  }

  if(r->order < 3) return;

  r->d3edrs3 = -6.0*ax/(r->rs[2]*r->rs[2]);
  if(p->nspin == XC_POLARIZED){
    if(ABS(r->zeta) == 1.0)
      d3fz = FLT_MAX;
    else
      d3fz = -4.0/27.0*(pow(1.0 + r->zeta,  -5.0/3.0) - pow(1.0 - r->zeta,  -5.0/3.0));

    r->d3edrs3 *= fz;
    r->d3edrs2z = 2.0*ax/(r->rs[1]*r->rs[2])*dfz;
    r->d3edrsz2 =    -ax/r->rs[2]           *d2fz;
    r->d3edz3   =     ax/r->rs[1]           *d3fz;
  }

  /* Relativistic corrections */
  /*  A. K. Rajagopal, J. Phys. C 11, L943 (1978).
      A. H. MacDonald and S. H. Vosko, J. Phys. C 12, 2977 (1979).
      E. Engel, S. Keller, A. Facco Bonetti, H. Müller, and R. M. Dreizler, Phys. Rev. A 52, 2750 (1995).
  */
  /*
  if(p->relativistic != 0){
    FLOAT beta, beta2, f1, f2, f3, phi;
  
    beta   = POW(3.0*M_PI*M_PI*dens, 1.0/3.0)/M_C;
    beta2  = beta*beta;
    f1     = sqrt(1.0 + beta2);
    f2     = asinh(beta);
    f3     = f1/beta - f2/beta2;
    phi    = 1.0 - 3.0/2.0*f3*f3;

    extmp *= phi;

    if(vrho != NULL){
      FLOAT dphidbeta, dbetadd;
      
      dphidbeta = 6.0/(beta2*beta2*beta)*
	(beta2 - beta*(2 + beta2)*f2/f1 + f2*f2);
      dbetadd = M_PI*M_PI/(POW(M_C, 3)*beta2);
    
      vrho[0]  = vrho[0]*phi + dens*extmp*dphidbeta*dbetadd;
    }
  }
  */
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_x) = {
  XC_LDA_X,
  XC_EXCHANGE,
  "Slater exchange",
  XC_FAMILY_LDA,
  "PAM Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930)\n"
  "F Bloch, Zeitschrift fuer Physik 57, 545 (1929)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC | XC_PROVIDES_KXC,
  lda_x_init,
  lda_x_end,
  work_lda
};

const XC(func_info_type) XC(func_info_lda_c_xalpha) = {
  XC_LDA_C_XALPHA,
  XC_CORRELATION,
  "Slater's Xalpha",
  XC_FAMILY_LDA,
  NULL,
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  lda_c_xalpha_init,
  lda_x_end,
  work_lda
};

