#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized 
 by 
   Perdew & Zunger
   Ortiz & Ballone
************************************************************************/

static func_type func_lda_c_pz = {
  XC_LDA_C_PZ,
  XC_CORRELATION,
  "Perdew & Zunger",
  "LDA",
  {"Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)", NULL}
};


void lda_c_pz_init(lda_type *p)
{
  p->func = &func_lda_c_pz;
}


static func_type func_lda_c_ob_pz = {
  XC_LDA_C_OB_PZ,
  XC_CORRELATION,
  "Ortiz & Ballone (PZ parametrization)",
  "LDA",
  {"Ortiz and Ballone, Phys. Rev. B 50, 1391 (1994)",
   "Ortiz and Ballone, Phys. Rev. B 56, 9970(E) (1997)",
   "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)", NULL}
};


void lda_c_ob_pz_init(lda_type *p)
{
  p->func = &func_lda_c_ob_pz;
}


/* Auxiliary functions to handle parametrizations */
static void ec_pot_low(int par, int i, double *rs, double *ec, double *pot)
{
  /* gamma[0], etc are PZ, and gamma[1] OB values */
  static double gamma[2][2] = {{-0.1423, -0.0843}, {-0.103756, -0.065951}};
  static double beta1[2][2] = {{ 1.0529,  1.3981}, { 0.56371,   1.11846}};
  static double beta2[2][2] = {{ 0.3334,  0.2611}, { 0.27358,   0.18797}};
  
  double b = 1.0/(1.0 + beta1[par][i]*rs[0] + beta2[par][i]*rs[1]);
  
  /* Eq. C3 */
  *ec  = gamma[par][i]*b;
  /* Eq. C4 */
  *pot = (*ec)*(1.0 + (7.0/6.0)*beta1[par][i]*rs[0] + 
		(4.0/3.0)*beta2[par][i]*rs[1])*b;
}


static void ec_pot_high(int par, int i, double *rs, double *ec, double *pot)
{
  /* the sign of c[1][0] and c[1][1] is diferent from [2], but is consistent
     with PWSCF. There is nothing in [3] about this, but I assume that PWSCF 
     is correct as it has the same sign as the PZ parametrizations */
  static double a[2][2] = {{ 0.0311,  0.015555}, { 0.031091,  0.015545}};
  static double b[2][2] = {{-0.048,  -0.0269},   {-0.046644, -0.025599}};
  static double c[2][2] = {{ 0.0020191519406228,  0.00069255121311694}, { 0.00419,  0.00329}};
  static double d[2][2] = {{-0.0116320663789130, -0.00480126353790614}, {-0.00983, -0.00300}};
  
  /* Eq. [1].C5 */
  *ec  = a[par][i]*rs[2] + b[par][i] + c[par][i]*rs[1]*rs[2] + d[par][i]*rs[1];
  /* Eq. [1].C6 */
  *pot = a[par][i]*rs[2] + (b[par][i] - a[par][i]/3.0) + (2.0/3.0)*c[par][i]*rs[1]*rs[2] + 
    (2.0*d[par][i] - c[par][i])*rs[1]/3.0;
}


/* the functional */
void lda_c_pz(lda_type *p, double rs_, double zeta, double *ec, double *vc)
{
  double rs[3];
  int func = p->func->number - XC_LDA_C_PZ;
  
  assert(func==0 || func==1);
  
  /* Wigner radius */
  rs[1] = rs_;
  rs[0] = sqrt(rs[1]);
  rs[2] = log(rs[1]);
  
  if(rs[1] >= 1.0)
    ec_pot_low (func, 0, rs, ec, &(vc[0]));
  else
    ec_pot_high(func, 0, rs, ec, &(vc[0]));
  
  if(p->nspin == XC_POLARIZED){
    double fz, fzp, ecp, vcp, x;
    
    fz  =  FZETA(zeta);
    fzp = DFZETA(zeta);
    
    if(rs[1] >= 1.0)
      ec_pot_low (func, 1, rs, &ecp, &vcp);
    else
      ec_pot_high(func, 1, rs, &ecp, &vcp);
    
    x = vc[0] + fz*(vcp - vc[0]) - zeta*(ecp - (*ec))*fzp;
    vc[0] = x + (ecp - (*ec))*fzp;
    vc[1] = x - (ecp - (*ec))*fzp;
    
    *ec += fz*(ecp - (*ec));
  }
}
