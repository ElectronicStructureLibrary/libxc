#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Wigner's parametrization from the low density limit

   [1] E.P. Wigner, Trans. Faraday Soc. 34, 678 (1938)
************************************************************************/

void lda_c_wigner(lda_type *p, double rs, double *ec, double *vc)
{
	static double a = -0.44, b = 7.8;
			
	double t = b + rs;
	
	*ec   =  a/t;
	vc[0] = -a/(t*t);                         /* now contains d ec/d rs */
	
	vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
	if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}


/************************************************************************
 Random Phase Approximation (RPA)

   [1] M. Gell-Mann & K.A. Brueckner, Phys. Rev. 106, 364 (1957)
************************************************************************/

void lda_c_rpa(lda_type *p, double rs, double *ec, double *vc)
{
	static double a = 0.0311, b = -0.047, c = 0.009, d = -0.017;

	double lrs = log(rs);

	*ec   = a*lrs + b + c*rs*lrs + d*rs;
	vc[0] = a/rs + c*(lrs + 1.0) + d;         /* now contains d ec/d rs */

	vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
	if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}


/************************************************************************
   [1] L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
     (can someone get me this paper, so I can find all coefficients?)
   [2] O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
************************************************************************/

static void hl_f(int func, int i, double rs, double *ec, double *vc)
{
	static const 
		double r[2][2] = {{21.0,   21.0},     /* HL unpolarized only*/
											{11.4,   15.9}};    /* GL */
	static const 
		double c[2][2] = {{ 0.0225, 0.0225},  /* HL unpolarized only */
											{ 0.0333, 0.0203}}; /* GL */

	double a, x, x2, x3;

	x   = rs/r[func][i];
	x2  = x*x;
  x3  = x2*x;

	a   = log(1.0 + 1.0/x);
	*ec = -c[func][i]*((1.0 + x3)*a - x2 + 0.5*x - 1.0/3.0);

	*vc = -c[func][i]*a;
}


void lda_c_hl(lda_type *p, double rs, double zeta, double *ec, double *vc)
{

	double ecp, vcp;
	int func = p->functional - XC_LDA_C_HL;

	/* sanity check */
	assert(func==0 || func==1);

	hl_f(func, 0, rs, &ecp, &vcp);

	if(p->nspin==XC_UNPOLARIZED){
		*ec   = ecp;
		vc[0] = vcp;

	}else{ /* XC_POLARIZED */
		double ecf, vcf, fz, dfz, t1;
		
		fz  =  FZETA(zeta);
		dfz = DFZETA(zeta);
		
		hl_f(func, 1, rs, &ecf, &vcf);
		
		*ec = ecp + (ecf - ecp)*fz;                  /* the energy    */
		t1  = vcp + (vcf - vcp)*fz;
		
		vc[0] = t1 + (ecf - ecp)*dfz*( 1.0 - zeta);  /* the potential */
		vc[1] = t1 + (ecf - ecp)*dfz*(-1.0 - zeta);
	}
}


/************************************************************************
 Slater's Xalpha functional

    Exc = alpha Ex
************************************************************************/

void lda_c_xalpha_init(lda_type *p, int nspin, int dim, int rel, double alpha)
{
	p->alpha = alpha;
	lda_x_init(p, nspin, dim, rel);
	p->functional = XC_LDA_C_XALPHA;
}

void lda_c_xalpha(lda_type *p, double *rho, double *ec, double *vc)
{
	int i;

	lda_x(p, rho, ec, vc);
	(*ec) *= p->alpha;
	for(i=0; i<p->nspin; i++) vc[i] *= p->alpha;
}


/************************************************************************
 LDA part of LYP's functional (translated from PWSCF)

   [1] C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
************************************************************************/

void lda_c_lyp(lda_type *p, double rs, double *ec, double *vc)
{
	static const double a=0.04918, b=0.132*2.87123400018819108;
  /* pi43=(4pi/3)^(1/3)     0.2533*pi43             0.349*pi43 */
	static const double c=0.408317561952371851, d=0.56258519195174803;
	
  double ecrs, ox;

  ecrs  = b*exp(-c*rs);
  ox    = 1.0/(1.0 + d*rs);
  *ec   = -a*ox*(1.0 + ecrs);
  vc[0] = (*ec) - rs/3.0*a*ox*(d*ox + ecrs*(d*ox + c));

	if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}

