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
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)"
};


void lda_c_pz_init(lda_type *p)
{
  p->func = &func_lda_c_pz;
}


static func_type func_lda_c_pz_mod = {
  XC_LDA_C_PZ_MOD,
  XC_CORRELATION,
  "Perdew & Zunger (Modified)",
  "LDA",
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)\n"
  "Modified to improve the matching between the low and high rs parts"
};

void lda_c_pz_mod_init(lda_type *p)
{
  p->func = &func_lda_c_pz_mod;
}


static func_type func_lda_c_ob_pz = {
  XC_LDA_C_OB_PZ,
  XC_CORRELATION,
  "Ortiz & Ballone (PZ parametrization)",
  "LDA",
  "Ortiz and Ballone, Phys. Rev. B 50, 1391 (1994)\n"
  "Ortiz and Ballone, Phys. Rev. B 56, 9970(E) (1997)\n"
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)"
};


void lda_c_ob_pz_init(lda_type *p)
{
  p->func = &func_lda_c_ob_pz;
}

typedef struct {
  double gamma[2];
  double beta1[2];
  double beta2[2];
  double a[2], b[2], c[2], d[2];
} pz_consts_type;

static pz_consts_type pz_consts[3] = {
  {    /* PZ Original */
    {-0.1423, -0.0843},  /* gamma */
    { 1.0529,  1.3981},  /* beta1 */
    { 0.3334,  0.2611},  /* beta2 */
    { 0.0311,  0.01555}, /*  a    */
    {-0.048,  -0.0269},  /*  b    */
    { 0.0020,  0.0007},  /*  c    */
    {-0.0116, -0.0048}   /*  d    */
  }, { /* PZ Modified */
    {-0.1423, -0.0843},   
    { 1.0529,  1.3981}, 
    { 0.3334,  0.2611}, 
    { 0.0311,  0.01555},
    {-0.048,  -0.0269},   
    { 0.0020191519406228,  0.00069255121311694},
    {-0.0116320663789130, -0.00480126353790614}
  }, { /* OB */
    {-0.103756, -0.065951},
    { 0.56371,   1.11846},
    { 0.27358,   0.18797},
    { 0.031091,  0.015545},
    {-0.046644, -0.025599},
    { 0.00419,   0.00329},  /* the sign of c[0] and c[1] is diferent from [2], but is consistent
			       with PWSCF. There is nothing in [3] about this, but I assume that PWSCF 
			       is correct as it has the same sign as the PZ parametrizations */
    {-0.00983,  -0.00300}
  }
};


/* Auxiliary functions to handle parametrizations */
static void ec_pot_low(pz_consts_type *X, int i, double *rs, double *ec, double *pot)
{
  double b = 1.0/(1.0 + X->beta1[i]*rs[0] + X->beta2[i]*rs[1]);
  
  /* Eq. C3 */
  *ec  = X->gamma[i]*b;

  /* Eq. C4 */
  *pot = (*ec)*(1.0 + (7.0/6.0)*X->beta1[i]*rs[0] + 
		(4.0/3.0)*X->beta2[i]*rs[1])*b;
}


static void ec_pot_high(pz_consts_type *X, int i, double *rs, double *ec, double *pot)
{
  /* Eq. [1].C5 */
  *ec  = X->a[i]*rs[2] + X->b[i] + X->c[i]*rs[1]*rs[2] + X->d[i]*rs[1];

  /* Eq. [1].C6 */
  *pot = X->a[i]*rs[2] + (X->b[i] - X->a[i]/3.0) + (2.0/3.0)*X->c[i]*rs[1]*rs[2] + 
    (2.0*X->d[i] - X->c[i])*rs[1]/3.0;
}


/* the functional */
void lda_c_pz(lda_type *p, double rs_, double zeta, double *ec, double *vc)
{
  double rs[3];
  int func = p->func->number - XC_LDA_C_PZ;
  
  assert(func==0 || func==1 || func==2);
  
  /* Wigner radius */
  rs[1] = rs_;
  rs[0] = sqrt(rs[1]);
  rs[2] = log(rs[1]);
  
  if(rs[1] >= 1.0)
    ec_pot_low (&pz_consts[func], 0, rs, ec, &(vc[0]));
  else
    ec_pot_high(&pz_consts[func], 0, rs, ec, &(vc[0]));
  
  if(p->nspin == XC_POLARIZED){
    double fz, fzp, ecp, vcp, x;
    
    fz  =  FZETA(zeta);
    fzp = DFZETA(zeta);
    
    if(rs[1] >= 1.0)
      ec_pot_low (&pz_consts[func], 1, rs, &ecp, &vcp);
    else
      ec_pot_high(&pz_consts[func], 1, rs, &ecp, &vcp);
    
    x = vc[0] + fz*(vcp - vc[0]) - zeta*(ecp - (*ec))*fzp;
    vc[0] = x + (ecp - (*ec))*fzp;
    vc[1] = x - (ecp - (*ec))*fzp;
    
    *ec += fz*(ecp - (*ec));
  }
}
