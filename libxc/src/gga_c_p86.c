#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Perdew 86 Generalized Gradient Approximation
 correlation functional.
************************************************************************/

void gga_c_p86_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_init(p->lda_aux, XC_LDA_C_PZ, p->nspin);
}

void gga_c_p86_end(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  free(p->lda_aux);
}

void gga_c_p86(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  double dens, zeta, dzdd[2], gdmt, ecunif, vcunif[2];
  double rs, DD, dDDdzeta, CC, CCinf, dCCdd;
  double Phi, dPhidd, dPhidgdmt;

  xc_lda(p->lda_aux, rho, &ecunif, vcunif, NULL);

  rho2dzeta(p->nspin, rho, &dens, &zeta);
  dzdd[0] =  (1.0 - zeta)/dens;
  dzdd[1] = -(1.0 + zeta)/dens;
    
  rs = RS(dens);

  /* get gdmt = |nabla n| */
  gdmt = sigma[0];
  if(p->nspin == XC_POLARIZED) gdmt += 2.0*sigma[1] + sigma[2];
  gdmt = sqrt(gdmt);


  { /* Equation [1].(4) */ 
    DD       = sqrt(pow(1.0 + zeta, 5.0/3.0) + pow(1.0 - zeta, 5.0/3.0))/M_SQRT2;
    dDDdzeta = 5.0/(3.0*4.0*DD)*(pow(1.0 + zeta, 2.0/3.0) - pow(1.0 - zeta, 2.0/3.0));
  }

  { /* Equation (6) of [1] */
    static const double alpha = 0.023266, beta = 7.389e-6, gamma = 8.723, delta = 0.472;
    static const double aa = 0.001667, bb = 0.002568;

    double rs2 = rs*rs, f1, f2, df1, df2, drsdd;

    f1    = bb + alpha*rs + beta*rs2;
    f2    = 1.0 + gamma*rs + delta*rs2 + 1.0e4*beta*rs*rs2;
    CC    = aa + f1/f2;
    CCinf = aa + bb;

    df1   = alpha + 2.0*beta*rs;
    df2   = gamma + 2.0*delta*rs + 3.0e4*beta*rs2;
    drsdd = -rs/(3.0*dens);
    dCCdd = (df1*f2 - f1*df2)/(f2*f2)*drsdd;
  }

  { /* Equation (9) of [1] */
    static const double ftilde = 1.745*0.11;

    double f1, f2, df1, df2;

    f1  = ftilde*(CCinf/CC);
    f2  = pow(dens, -7.0/6.0);
    Phi = f1*gdmt*f2;

    df1 = -f1/(CC)*dCCdd;
    df2 = -7.0/6.0*pow(dens, -13.0/6.0);
    dPhidd    = gdmt*(df1*f2 + f1*df2);
    dPhidgdmt = f1*f2;
  }

  { /* Equation [1].(8) */
    double gdmt2;
    double f1, f2, f3, df1, df1dgdmt, df2, df3, df3dgdmt;

    gdmt2 = gdmt*gdmt;

    f1 = exp(-Phi);
    f2 = pow(dens, -4.0/3.0);
    f3 = f1*CC*gdmt2*f2;

    df1      = -f1*dPhidd;
    df1dgdmt = -f1*dPhidgdmt;
    df2      = -4.0/3.0*pow(dens, -7.0/3.0);
    df3      = gdmt2*(df1*CC*f2 + f1*dCCdd*f2 + f1*CC*df2);
    df3dgdmt = CC*f2*(df1dgdmt*gdmt2 + f1*2.0*gdmt);

    *e = ecunif + f3/(DD*dens);

    vrho[0]   = vcunif[0] + (df3 - (f3/DD)*dDDdzeta*dzdd[0])/DD;
    vsigma[0] = df3dgdmt/(DD*2.0*gdmt);

    if(p->nspin == XC_POLARIZED){
      vrho[1]   = vcunif[1] + (df3 - (f3/DD)*dDDdzeta*dzdd[1])/DD;
      vsigma[1] = 2.0*vsigma[0];
      vsigma[2] =     vsigma[0];
    }
  }
}

const xc_func_info_type func_info_gga_c_p86 = {
  XC_GGA_C_P86,
  XC_CORRELATION,
  "Perdew 86",
  XC_FAMILY_GGA,
  "J.P.Perdew, Phys. Rev. B 33, 8822 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_c_p86_init,
  gga_c_p86_end,   /* we can use the same as exchange here */
  NULL,            /* this is not an LDA                   */
  gga_c_p86
};
