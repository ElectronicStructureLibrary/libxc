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
  {"J.P. Perdew and Y.Wang, Phys. Rev. B 45, 13244 (1992)", NULL}
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
  {"Ortiz and Ballone, Phys. Rev. B 50, 1391 (1994)",
   "Ortiz and Ballone, Phys. Rev. B 56, 9970(E) (1997)",
   "J.P. Perdew and Y.Wang, Phys. Rev. B 45, 13244 (1992)", NULL}
};


void lda_c_ob_pw_init(lda_type *p)
{
  p->func = &func_lda_c_ob_pw;
}


/* Function g defined by Eq. 10 of the original paper,
   and it's derivative with respect to rs, Eq. A5 */
static void g(int func, int k, double *rs, double *f, double *dfdrs)
{
  static double a[2][3]     = {{0.031091, 0.015545, 0.016887},   /* PZ */
			       {0.031091, 0.015545, 0.016887}};  /* OB */
  static double alpha[2][3] = {{0.21370,  0.20548,  0.11125},    /* PZ */
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
  
  double q0, q1, q1p, b;
  
  b = beta[func][k][0]*rs[0] + beta[func][k][1]*rs[1] + 
    beta[func][k][2]*rs[0]*rs[1] + beta[func][k][3]*rs[2];
  
  /* the function */
  *f = -2.0*a[func][k]*(1.0 + alpha[func][k]*rs[1])*log(1.0 + 1.0/(2.0*a[func][k]*b));
  
  /* and now the derivative */
  q0     = -2.0*a[func][k]*(1.0 + alpha[func][k] * rs[1]);
  q1     =  2.0*a[func][k]*b;
  q1p    = a[func][k]*(beta[func][k][0]/rs[0] + 2.0*beta[func][k][1] + 
		       3.0*beta[func][k][2]*rs[0] + 4.0*beta[func][k][3]*rs[1]);
  *dfdrs = -2.0*a[func][k]*alpha[func][k]*log(1.0 + 1.0/q1) - (q0*q1p)/(q1*q1 + q1);
}


/* the functional */
void lda_c_pw(lda_type *p, double rs_, double zeta, double *ec, double *vc)
{
  double rs[3], decdrs;
  int func = p->func->number - XC_LDA_C_PW;
  
  assert(func==0 || func==1);
  
  /* Wigner radius */
  rs[1] = rs_;
  rs[0] = sqrt(rs[1]);
  rs[2] = rs[1]*rs[1];
  
  g(func, 0, rs, ec, &decdrs);
  
  if(p->nspin == XC_UNPOLARIZED){
    vc[0] = (*ec) - (rs[1]/3.0)*decdrs;
    
  }else{
    static double fz20 = 1.709921;          /* = (d^2f/dz^2)(z=0) */
    
    double fz, fpz, z4, ec1, alphac, dec1drs, dalphacdrs, decdz;
    
    fz  =  FZETA(zeta);
    fpz = DFZETA(zeta);
    z4  = pow(zeta, 4);
    
    g(func, 1, rs, &ec1,    &dec1drs);
    g(func, 2, rs, &alphac, &dalphacdrs);
    
    decdrs = decdrs*(1.0 - fz*z4) + dec1drs*fz*z4 + dalphacdrs*fz*(1.0 - z4)/fz20;
    decdz  = 4.0*pow(zeta, 3)*fz*(ec1 - (*ec) - alphac/fz20) + 
      fpz*(z4*(ec1 - (*ec)) + (1.0 - z4)*alphac/fz20);
    
    *ec += alphac*fz*(1.0 - z4)/fz20 + (ec1 - (*ec))*fz*z4;
    
    vc[0] = (*ec) - (rs[1]/3.0)*decdrs - (zeta - 1.0)*decdz;
    vc[1] = (*ec) - (rs[1]/3.0)*decdrs - (zeta + 1.0)*decdz;
  }
}
