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

 If irel is not zero, relativistic correction factors have to be applied.
 These are however not implemented in 2D, so nothing is done in that case.

 WARNING: Check that the relativistic corrections are OK for the potential
          in the spin polarized case.
************************************************************************/

static func_type func_lda_x = {
  XC_LDA_X,
  XC_EXCHANGE,
  "Slater exchange",
  "LDA",
  NULL
};


void lda_x_init(lda_type *p, int nspin, int dim)
{
  p->func = &func_lda_x;

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  assert(dim>=2 && dim<=3);
  p->dim = dim;
}


void lda_x(lda_type *p, double *rho, double *ex, double *vx, double *fx)
{
  static double a_x[3] = {-1.0, -1.06384608107049, -0.738558766382022};

  double dens, alpha, factor;
  int i;
  
  assert(p!=NULL);
  
  dens = 0.0;
  for(i=0; i<p->nspin; i++) dens += rho[i];

  if(dens <= MIN_DENS){
    *ex = 0.0;
    for(i=0; i<p->nspin; i++) vx[i] = 0.0;
    if(fx != NULL)
      for(i=0; i<p->nspin*p->nspin; i++) fx[i] = 0.0;
    return;
  }
  
  alpha = (p->dim + 1.0)/p->dim;
  factor = (p->nspin == XC_UNPOLARIZED) ? 1.0 : pow(2.0, alpha)/2.0;
  factor *= a_x[p->dim-1];

  *ex = 0.0;
  if(fx != NULL)
    for(i=0; i<p->nspin*p->nspin; i++) fx[i] = 0.0;

  for(i=0; i<p->nspin; i++){
    *ex   += factor*pow(rho[i], alpha);
    vx[i]  = factor*alpha*pow(rho[i], alpha - 1.0);
    if(fx!=NULL && rho[i]>0)
      fx __(i,i) = factor*alpha*(alpha - 1.0)*pow(rho[i], alpha - 2.0);
  }
  *ex /= dens;
}
