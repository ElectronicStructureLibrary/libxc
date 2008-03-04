/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized 
 by 
   Perdew & Zunger
   Ortiz & Ballone
************************************************************************/

#define XC_LDA_C_PZ       9   /* Perdew & Zunger              */
#define XC_LDA_C_PZ_MOD  10   /* Perdew & Zunger (Modified)   */
#define XC_LDA_C_OB_PZ   11   /* Ortiz & Ballone (PZ)         */

typedef struct {
  FLOAT gamma[2];
  FLOAT beta1[2];
  FLOAT beta2[2];
  FLOAT a[2], b[2], c[2], d[2];
} pz_consts_type;

static pz_consts_type
pz_consts[3] = {
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
static void
ec_pot_low(pz_consts_type *X, int i, FLOAT *rs, FLOAT *ec, FLOAT *pot)
{
  FLOAT b = 1.0/(1.0 + X->beta1[i]*rs[0] + X->beta2[i]*rs[1]);
  
  /* Eq. C3 */
  *ec  = X->gamma[i]*b;

  /* Eq. C4 */
  *pot = (*ec)*(1.0 + (7.0/6.0)*X->beta1[i]*rs[0] + 
		(4.0/3.0)*X->beta2[i]*rs[1])*b;
}


static void 
ec_pot_high(pz_consts_type *X, int i, FLOAT *rs, FLOAT *ec, FLOAT *pot)
{
  /* Eq. [1].C5 */
  *ec  = X->a[i]*rs[2] + X->b[i] + X->c[i]*rs[1]*rs[2] + X->d[i]*rs[1];

  /* Eq. [1].C6 */
  *pot = X->a[i]*rs[2] + (X->b[i] - X->a[i]/3.0) + (2.0/3.0)*X->c[i]*rs[1]*rs[2] + 
    (2.0*X->d[i] - X->c[i])*rs[1]/3.0;
}


/* the functional */
static void
lda_c_pz(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  FLOAT dens, zeta, rs[3];
  int func = p->info->number - XC_LDA_C_PZ;
  
  assert(func==0 || func==1 || func==2);
  
  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);

  /* Wigner radius */
  rs[1] = RS(dens);
  rs[0] = sqrt(rs[1]);
  rs[2] = log(rs[1]);
  
  if(rs[1] >= 1.0)
    ec_pot_low (&pz_consts[func], 0, rs, ec, &(vc[0]));
  else
    ec_pot_high(&pz_consts[func], 0, rs, ec, &(vc[0]));
  
  if(p->nspin == XC_POLARIZED){
    FLOAT fz, fzp, ecp, vcp, x;
    
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

const XC(func_info_type) XC(func_info_lda_c_pz) = {
  XC_LDA_C_PZ,
  XC_CORRELATION,
  "Perdew & Zunger",
  XC_FAMILY_LDA,
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,
  NULL,
  lda_c_pz
};

const XC(func_info_type) XC(func_info_lda_c_pz_mod) = {
  XC_LDA_C_PZ_MOD,
  XC_CORRELATION,
  "Perdew & Zunger (Modified)",
  XC_FAMILY_LDA,
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)\n"
  "Modified to improve the matching between the low and high rs parts",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,
  NULL,
  lda_c_pz
};

const XC(func_info_type) XC(func_info_lda_c_ob_pz) = {
  XC_LDA_C_OB_PZ,
  XC_CORRELATION,
  "Ortiz & Ballone (PZ parametrization)",
  XC_FAMILY_LDA,
  "G Ortiz and P Ballone, Phys. Rev. B 50, 1391 (1994)\n"
  "G Ortiz and P Ballone, Phys. Rev. B 56, 9970(E) (1997)\n"
  "Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,
  NULL,
  lda_c_pz
};
