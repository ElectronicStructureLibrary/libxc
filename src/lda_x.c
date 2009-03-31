/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

#define XC_LDA_X  1   /* Exchange                     */

static inline void 
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  const FLOAT ax = -0.458165293283142893475554485052; /* -3/4*POW(3/(2*M_PI), 2/3) */
  FLOAT fz, dfz, d2fz, d3fz;

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
  NULL, NULL,
  work_lda
};


void XC(lda_x_init)(XC(lda_type) *p, int nspin, int dim, int irel)
{
  p->info = &XC(func_info_lda_x);

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  assert(dim>=2 && dim<=3);
  assert(irel == 0 || (dim==3 && nspin==XC_UNPOLARIZED));

  p->dim = dim;
  p->relativistic = irel;
  p->nspin = nspin;
}


