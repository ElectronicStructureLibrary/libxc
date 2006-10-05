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

void lda_x(void *p_, double *rho, double *ex, double *vx, double *fx)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a_x[3] = {-1.0, -1.06384608107049, -0.738558766382022};
  double dens, alpha, factor, beta, phi, DphiDdens;
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
  
  alpha   = (p->dim + 1.0)/p->dim;
  factor  = (p->nspin == XC_UNPOLARIZED) ? 1.0 : pow(2.0, alpha)/2.0;
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

  if(p->relativistic != 0){
    /* Relativistic corrections */
    beta = pow(3.0*pow(M_PI, 2.0)*dens, 1.0/3.0)/M_C;
    phi = 1.0 - 3.0/2.0*pow(sqrt(1 + pow(beta, 2.0))/beta - asinh(beta)/pow(beta, 2.0), 2.0);
    DphiDdens = -2.0/dens/pow(beta, 2.0) * 
      (-1.0 + asinh(beta)*(pow(beta, 2.0) + 2)/(beta*sqrt(1.0 + pow(beta, 2.0))) -
       pow(asinh(beta), 2.0)/pow(beta, 2.0));

    vx[0] = vx[0]*phi + dens*DphiDdens*(*ex);
    *ex *= phi;
  };
}


const xc_func_info_type func_info_lda_x = {
  XC_LDA_X,
  XC_EXCHANGE,
  "Slater exchange",
  XC_FAMILY_LDA,
  "P.A.M. Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930)\n"
  "F. Bloch, Zeitschrift für Physik 57, 545 (1929)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,
  NULL,
  lda_x
};


void xc_lda_x_init(xc_lda_type *p, int nspin, int dim, int irel)
{
  p->info = &func_info_lda_x;

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  assert(dim>=2 && dim<=3);
  assert(irel == 0 || (dim==3 && nspin==XC_UNPOLARIZED));

  p->dim = dim;
  p->relativistic = irel;
  p->nspin = nspin;
}


