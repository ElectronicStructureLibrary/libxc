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

/************************************************************************
 Calculates the exchange energy and exchange potential for a homogeneous
 electron gas (LDA for exchange).
 It works in both the 3D and 2D cases
 (1D not yet implemented, although one) only has to find out the a_x constant,
 but I do not know what is the value of \int_0^\infty sin**2(x)/x**3 )

 The basic formulae are (Hartree atomic units are assumed):

    p = ((dim+1)/dim)
    ex(n) = a_x*n**(1/dim)
    ex(n,z) = ex(n)*f(z)
    vx_up(n, z) = ex(n)*( p*f(z) + (df/dz)(z)*(1-z) )
    vx_do(n, z) = ex(n)*( p*f(z) - (df/dz)(z)*(1+z) )
    f(z) = (1/2)*( (1+z)**p + (1-z)**p)
    a_x = -(3/(4*pi))*(3*pi**2)**(1/3) in 3D
    a_x = -(4/3)*sqrt(2/pi) in 2D
    a_x = -(1/2) * \int_0^\infty (sin(x))**2/x**3

 If irel is not zero, a relativistic correction factor is applied.
 This factor can only be aplied in 3D and for the spin-unpolarized case.
************************************************************************/

#define XC_LDA_X  1   /* Exchange                     */

static void lda_x(const void *p_, const FLOAT *rho, FLOAT *ex, FLOAT *vx, FLOAT *fx)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  static FLOAT a_x[3] = {-1.0, -1.06384608107049, -0.738558766382022};
  FLOAT dens, extmp, alpha, factor;
  int i;
  
  assert(p!=NULL);
  
  dens = 0.0;
  for(i=0; i<p->nspin; i++) dens += rho[i];

  alpha   = (p->dim + 1.0)/p->dim;
  factor  = (p->nspin == XC_UNPOLARIZED) ? 1.0 : POW(2.0, alpha)/2.0;
  factor *= a_x[p->dim-1];

  extmp = 0.0;
  for(i=0; i<p->nspin; i++){
    extmp += factor*POW(rho[i], alpha)/dens;

    if(vx != NULL)
      vx[i]  = factor*alpha*POW(rho[i], alpha - 1.0);

    if(fx!=NULL && rho[i]>0){
      int js = (i==0) ? 0 : 2;
      fx[js] = factor*alpha*(alpha - 1.0)*POW(rho[i], alpha - 2.0);
    }
  }

  /* Relativistic corrections */
  if(p->relativistic != 0){
    FLOAT beta, beta2, f1, f2, f3, phi;
  
    beta   = POW(3.0*M_PI*M_PI*dens, 1.0/3.0)/M_C;
    beta2  = beta*beta;
    f1     = sqrt(1.0 + beta2);
    f2     = asinh(beta);
    f3     = f1/beta - f2/beta2;
    phi    = 1.0 - 3.0/2.0*f3*f3;

    extmp *= phi;

    if(vx != NULL){
      FLOAT dphidbeta, dbetadd;
      
      dphidbeta = 6.0/(beta2*beta2*beta)*
	(beta2 - beta*(2 + beta2)*f2/f1 + f2*f2);
      dbetadd = M_PI*M_PI/(POW(M_C, 3)*beta2);
    
      vx[0]  = vx[0]*phi + dens*extmp*dphidbeta*dbetadd;
    }

    /* WARNING - missing fxc */
  }

  if(ex != NULL)
    *ex = extmp;
}


const XC(func_info_type) XC(func_info_lda_x) = {
  XC_LDA_X,
  XC_EXCHANGE,
  "Slater exchange",
  XC_FAMILY_LDA,
  "PAM Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930)\n"
  "F Bloch, Zeitschrift fuer Physik 57, 545 (1929)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,
  NULL,
  lda_x
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


