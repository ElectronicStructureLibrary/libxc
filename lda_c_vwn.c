#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 LDA parametrization of Vosko, Wilk & Nusair
************************************************************************/

static func_type func_lda_c_vwn = {
  XC_LDA_C_VWN,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair",
  "LDA",
  "S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)"
};

/* some constants */
static double  A[3] = { 0.0621841, 0.0310907, -0.033774 }; /* CHECK */
static double  b[3] = { 3.72744,   7.06042,    1.131071 };
static double  c[3] = {12.9352,   18.0578,    13.0045   };
static double x0[3] = {-0.10498,  -0.32500,   -0.0047584};
static double  Q[3] = { 0.0,       0.0,        0.0      };
static double fpp   = 0.0;

/* initialization */
void lda_c_vwn_init(lda_type *p)
{
  int i;
  
  p->func = &func_lda_c_vwn;

  for(i=0; i<3; i++){
    Q[i] = sqrt(4.0*c[i] - b[i]*b[i]);
    A[i] = A[i]/2.0; /* I believe these numbers are in Rydberg, so I convert to Hartree*/
  }
  fpp = 4.0/(9.0*(pow(2.0, 1.0/3.0)-1));
}

/* useful functions */
void ec_i(int i, double x, double *ec, double *decdrs)
{
  double f1, f2, f3, fx, qx, xx0, tx, tt;
  
  f1  = 2.0*b[i]/Q[i];
  f2  = b[i]*x0[i]/(x0[i]*x0[i] + b[i]*x0[i] + c[i]);
  f3  = 2.0*(2.0*x0[i] + b[i])/Q[i];
  fx  = x*x + b[i]*x + c[i];
  qx  = atan(Q[i]/(2.0*x + b[i]));
  xx0 = x - x0[i];
  
  *ec = A[i]*(log(x*x/fx) + f1*qx - f2*(log(xx0*xx0/fx) + f3*qx));
  
  tx  = 2.0*x + b[i];
  tt  = tx*tx + Q[i]*Q[i];
  *decdrs = A[i]*(2.0/x - tx/fx - 4.0*b[i]/tt -
		  f2*(2.0/xx0 - tx/fx - 4.0*(2.0*x0[i] + b[i])/tt));
}

/* the functional */
void lda_c_vwn(lda_type *p, double rs_, double zeta, double *ec, double *vc)
{
	double rs[2], ec_1, dec_1;

	assert(p!=NULL);

	/* Wigner radius */
	rs[1] = rs_;          /* rs          */
	rs[0] = sqrt(rs[1]);  /* sqrt(rs)    */

	ec_i(0, rs[0], &ec_1, &dec_1);

	if(p->nspin == XC_UNPOLARIZED){
		*ec   = ec_1;
		vc[0] = ec_1 - rs[0]/6.0*dec_1;
	}else{
		double ec_2, ec_3, dec_2, dec_3, fz, dfz, decdx, decdz;
		double t1, t2, z3, z4;

		ec_i(1, rs[0], &ec_2, &dec_2);
		ec_i(2, rs[0], &ec_3, &dec_3);

		fz  =  FZETA(zeta);
		dfz = DFZETA(zeta);
		
		z3 = pow(zeta, 3);
		z4 = z3*zeta;
		t1 = (fz/fpp)*(1.0 - z4);
		t2 = fz*z4;

		*ec   =  ec_1 +  ec_3*t1 + ( ec_2 -  ec_1)*t2;
		decdx = dec_1 + dec_3*t1 + (dec_2 - dec_1)*t2;
		decdz = (ec_3/fpp)*(dfz*(1.0 - z4) - 4.0*fz*z3) +
			(ec_2 - ec_1)*(dfz*z4 + 4.0*fz*z3);

		t1 = *ec - rs[0]/6.0*decdx;

		vc[0] = t1 + (1.0 - zeta)*decdz;
		vc[1] = t1 - (1.0 + zeta)*decdz;
	}
	
}
