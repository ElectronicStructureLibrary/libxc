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

void lda_x_init(lda_type *p, int nspin, int dim, int rel)
{
	p->functional = XC_LDA_X;

	assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
	p->nspin = nspin;

	assert(rel==XC_NON_RELATIVISTIC || rel==XC_RELATIVISTIC);
	p->relativistic = rel;

	assert(dim>=2 && dim<=3);
	p->dim = dim;
}


void lda_x(lda_type *p, double *rho, double *ex, double *vx)
{
	static double a_x[3] = {-1.0, -1.06384608107049, -0.738558766382022};

	double dens, zeta, alpha;
	int i;

	assert(p!=NULL);

	rho2dzeta(p->nspin, rho, &dens, &zeta);

	if(dens <= MIN_DENS){
		*ex = 0.0;
		for(i=0; i<p->nspin; i++) vx[i] = 0.0;
		return;
	}

	alpha = (p->dim + 1.0)/p->dim;
	*ex   = a_x[p->dim-1]*pow(dens, 1.0/p->dim);

	if(p->nspin == XC_UNPOLARIZED)
		vx[0] = alpha*(*ex);
	else{ /* XC_POLARIZED */
		double fz, dfdz;

		fz   = 0.5*(pow(1+zeta, alpha) + pow(1-zeta, alpha));
		dfdz = 0.5*alpha*(pow(1+zeta, alpha-1) - pow(1-zeta, alpha-1));

		vx[0] = (*ex)*(alpha*fz + dfdz*(1 - zeta));
		vx[1] = (*ex)*(alpha*fz - dfdz*(1 + zeta));
		(*ex) = (*ex)*fz; /* Now it is the correct polarized result. */
	}

	/* Relativistic corrections (I believe that they are rather useless). */
	if(p->relativistic == XC_RELATIVISTIC){
		double rs, beta, sb, alb, t;

		rs = RS(dens);
		beta = 0.014/rs;
		sb = sqrt(1 + beta*beta);
		alb = log(beta + sb);
		t = (beta*sb - alb)/(beta*beta);
		*ex = (*ex)*(1.0 - 1.5*t*t);
		for(i=0; i<p->nspin; i++) vx[i] = vx[i]*(-0.5 + 1.5*alb/(beta*sb));
	}

}
