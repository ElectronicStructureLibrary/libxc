#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized 
 by 
   J.P. Perdew & Y.Wang
   Ortiz & Ballone
************************************************************************/

static func_type func_lda_c_pw = {
  XC_LDA_C_PW,
  XC_CORRELATION,
  "Perdew & Wang",
  "LDA",
  "J.P. Perdew and Y.Wang, Phys. Rev. B 45, 13244 (1992)"
};


void lda_c_pw_init(lda_type *p)
{
  p->func = &func_lda_c_pw;
}


static func_type func_lda_c_ob_pw = {
  XC_LDA_C_OB_PW,
  XC_CORRELATION,
  "Ortiz & Ballone (PW parametrization)",
  "LDA",
  "Ortiz and Ballone, Phys. Rev. B 50, 1391 (1994)\n"
  "Ortiz and Ballone, Phys. Rev. B 56, 9970(E) (1997)\n"
  "J.P. Perdew and Y.Wang, Phys. Rev. B 45, 13244 (1992)"
};


void lda_c_ob_pw_init(lda_type *p)
{
  p->func = &func_lda_c_ob_pw;
}


/* Function g defined by Eq. 10 of the original paper,
   and it's derivative with respect to rs, Eq. A5 */
static void g(int func, int k, double *rs, double *f, double *dfdrs, double *d2fdrs2)
{
  static double a[2][3]     = {{0.031091, 0.015545, 0.016887},   /* PW */
			       {0.031091, 0.015545, 0.016887}};  /* OB */
  static double alpha[2][3] = {{0.21370,  0.20548,  0.11125},    /* PW */
			       {0.026481, 0.022465, 0.11125}};   /* OB */
  static double beta[2][3][4] = {
    {
      { 7.5957,  3.5876,   1.6382,  0.49294}, /* PZ */
      {14.1189,  6.1977,   3.3662,  0.62517},
      {10.357,   3.6231,   0.88026, 0.49671}
    },{
      { 7.5957,  3.5876,  -0.46647, 0.13354}, /* OB */
      {14.1189,  6.1977,  -0.56043, 0.11313},
      {10.357,   3.6231,   0.88026, 0.49671}
    }};
  
  double q0, q1, q1p, b, aux;
  
  b = beta[func][k][0]*rs[0] + beta[func][k][1]*rs[1] + 
    beta[func][k][2]*rs[0]*rs[1] + beta[func][k][3]*rs[2];
  
  q0     = -2.0*a[func][k]*(1.0 + alpha[func][k] * rs[1]);
  q1     =  2.0*a[func][k]*b;

  /* the function */
  *f = q0*log(1.0 + 1.0/q1);
  
  /* and now the derivative */
  aux = 1.0/(q1*q1 + q1);
  q1p = a[func][k]*(beta[func][k][0]/rs[0] + 2.0*beta[func][k][1] + 
		    3.0*beta[func][k][2]*rs[0] + 4.0*beta[func][k][3]*rs[1]);

  *dfdrs = -2.0*a[func][k]*alpha[func][k]*log(1.0 + 1.0/q1) - (q0*q1p)*aux;

  if(d2fdrs2 != NULL){
    double q1pp;

    q1pp = a[func][k]*(-beta[func][k][0]/(2.0*rs[0]*rs[1]) +
		       3.0*beta[func][k][2]/(2.0*rs[0]) + 4.0*beta[func][k][3]);

    *d2fdrs2 = aux*(4.0*a[func][k]*alpha[func][k]*q1p - q0*q1pp + q0*q1p*q1p*(2.0*q1 + 1.0)*aux);
  }
}


/* the functional */
void lda_c_pw(lda_type *p, double rs_, double dens, double zeta, double *ec, double *vc, double *fc)
{
  double rs[3], Dec_Drs, D2ec_Drs2, *dp;
  int func = p->func->number - XC_LDA_C_PW;
  
  assert(func==0 || func==1);
  
  /* Wigner radius */
  rs[1] = rs_;
  rs[0] = sqrt(rs[1]);
  rs[2] = rs[1]*rs[1];
  
  dp = (fc == NULL) ? NULL : (&D2ec_Drs2);
  g(func, 0, rs, ec, &Dec_Drs, dp);
  
  if(p->nspin == XC_UNPOLARIZED){
    vc[0] = (*ec) - (rs[1]/3.0)*Dec_Drs;
    
    if(fc != NULL){
      double Drs = -(4.0*M_PI/9.0)*rs[2]*rs[2];
      fc __(0,0) = (2.0*Dec_Drs - rs[1]*D2ec_Drs2)*Drs/3.0;
    }

  }else{
    static double fz20 = 1.709921;          /* = (d^2f/dz^2)(z=0) */
    
    double fz, fpz, z4, ec1, alphac;
    double ec0, Dec0_Drs, Dec1_Drs, D2ec1_Drs2, Dalphac_Drs, D2alphac_Drs2, Dec_Dz;
    
    fz  =  FZETA(zeta);
    fpz = DFZETA(zeta);
    z4  = pow(zeta, 4);

    dp = (fc == NULL) ? NULL : (&D2ec1_Drs2);
    g(func, 1, rs, &ec1, &Dec1_Drs, dp);

    dp = (fc == NULL) ? NULL : (&D2alphac_Drs2);
    g(func, 2, rs, &alphac, &Dalphac_Drs, dp);
    
    /* save copies that will be needed later */
    ec0      = (*ec);
    Dec0_Drs = Dec_Drs;  

    *ec     =  ec0 + z4*fz*(ec1 - ec0 - alphac/fz20) + fz*alphac/fz20;

    Dec_Drs = Dec0_Drs + z4*fz*(Dec1_Drs - Dec0_Drs - Dalphac_Drs/fz20) + fz*Dalphac_Drs/fz20;
    Dec_Dz  = 4.0*pow(zeta, 3)*fz*(ec1 - ec0 - alphac/fz20) + 
      fpz*z4*(ec1 - ec0 - alphac/fz20) + fpz*alphac/fz20;
    
    vc[0] = (*ec) - (rs[1]/3.0)*Dec_Drs - (zeta - 1.0)*Dec_Dz;
    vc[1] = (*ec) - (rs[1]/3.0)*Dec_Drs - (zeta + 1.0)*Dec_Dz;
    
    if(fc != NULL){
      double tmp, fppz, D2ec_Dz2, D2ec_DrsDz;
      double Drs = -(4.0*M_PI/9.0)*rs[2]*rs[2];
      int i;

      fppz = 0.0;
      if(zeta > -1.0) fppz += pow(1.0 + zeta, -2.0/3.0);
      if(zeta <  1.0) fppz += pow(1.0 - zeta, -2.0/3.0);
      fppz *= 4.0/(9.0*FZETAFACTOR);

      D2ec_Drs2  = D2ec_Drs2 + z4*fz*(D2ec1_Drs2 - D2ec_Drs2 - D2alphac_Drs2/fz20) + fz*D2alphac_Drs2/fz20;

      D2ec_DrsDz = 4.0*pow(zeta, 3)*fz*(Dec1_Drs - Dec0_Drs - Dalphac_Drs/fz20) + 
	fpz*(z4*(Dec1_Drs - Dec0_Drs) + (1.0 - z4)*Dalphac_Drs/fz20);

      D2ec_Dz2   = 4.0*pow(zeta, 2)*(ec1 - ec0 - alphac/fz20)*(3.0*fz + 2.0*zeta*fpz) + 
	fppz*(z4*(ec1 - ec0) + (1.0 - z4)*alphac/fz20);
    
      tmp = (2.0*Dec_Drs - rs[1]*D2ec_Drs2)*Drs/3.0;
      for(i=0; i<4; i++) fc[i]  = tmp;

      fc __(0,0) += -D2ec_DrsDz*(zeta - 1.0)*Drs 
	+ (zeta - 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta - 1.0)*D2ec_Dz2);

      fc __(1,1) += -D2ec_DrsDz*(zeta + 1.0)*Drs 
	+ (zeta + 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta + 1.0)*D2ec_Dz2);

      fc __(0,1) += -D2ec_DrsDz*(zeta - 1.0)*Drs 
	+ (zeta + 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta - 1.0)*D2ec_Dz2);

      fc __(1,0)  = fc __(0,1);
    }
  }
}
