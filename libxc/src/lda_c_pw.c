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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Correlation energy per-particle and potential of a HEG as parameterized 
 by 
   J.P. Perdew & Y. Wang
   Ortiz & Ballone

Note that the PW modified, corresponds to the version of PW used in the 
original PBE routine. This amounts to adding some more digits in some of
the constants of PW.
************************************************************************/

#define XC_LDA_C_PW     12   /* Perdew & Wang                */
#define XC_LDA_C_PW_MOD 13   /* Perdew & Wang (Modified)     */
#define XC_LDA_C_OB_PW  14   /* Ortiz & Ballone (PW)         */


/* Function g defined by Eq. 10 of the original paper,
   and it's derivative with respect to rs, Eq. A5 */
static void g(int func, int k, FLOAT *rs, FLOAT *f, FLOAT *dfdrs, FLOAT *d2fdrs2)
{
  static FLOAT a[3][3]     = 
    {
      {0.031091,  0.015545,   0.016887},    /* PW */
      {0.0310907, 0.01554535, 0.0168869},   /* PW (modified) */
      {0.031091,  0.015545,   0.016887}     /* OB */
  }; 
  static FLOAT alpha[3][3] = 
    {
      {0.21370,  0.20548,  0.11125},    /* PW */
      {0.21370,  0.20548,  0.11125},    /* PW (modified) */
      {0.026481, 0.022465, 0.11125}     /* OB */
    };
  static FLOAT beta[3][3][4] = {
    {
      { 7.5957,  3.5876,   1.6382,  0.49294}, /* PW */
      {14.1189,  6.1977,   3.3662,  0.62517},
      {10.357,   3.6231,   0.88026, 0.49671}
    },{
      { 7.5957,  3.5876,   1.6382,  0.49294}, /* PW (modified) */
      {14.1189,  6.1977,   3.3662,  0.62517},
      {10.357,   3.6231,   0.88026, 0.49671}
    },{
      { 7.5957,  3.5876,  -0.46647, 0.13354}, /* OB */
      {14.1189,  6.1977,  -0.56043, 0.11313},
      {10.357,   3.6231,   0.88026, 0.49671}
    }};
  
  FLOAT q0, q1, q1p, b, aux;
  
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
    FLOAT q1pp;

    q1pp = a[func][k]*(-beta[func][k][0]/(2.0*rs[0]*rs[1]) +
		       3.0*beta[func][k][2]/(2.0*rs[0]) + 4.0*beta[func][k][3]);

    *d2fdrs2 = aux*(4.0*a[func][k]*alpha[func][k]*q1p - q0*q1pp + q0*q1p*q1p*(2.0*q1 + 1.0)*aux);
  }
}


/* the functional */
void lda_c_pw(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  FLOAT dens, zeta;
  FLOAT rs[3], Dec_Drs, D2ec_Drs2, ec0, *dp;
  int func = p->info->number - XC_LDA_C_PW;
  
  assert(func==0 || func==1 || func==2);
  
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  /* Wigner radius */
  rs[1] = RS(dens);
  rs[0] = sqrt(rs[1]);
  rs[2] = rs[1]*rs[1];
  
  /* ec(rs, 0) */
  dp = (fc == NULL) ? NULL : (&D2ec_Drs2);
  g(func, 0, rs, &ec0, &Dec_Drs, dp);
  
  if(p->nspin == XC_UNPOLARIZED){
    if(ec != NULL)
      *ec = ec0;

    if(vc != NULL)
      vc[0] = (*ec) - (rs[1]/3.0)*Dec_Drs;
    
    if(fc != NULL){
      FLOAT Drs = -(4.0*M_PI/9.0)*rs[2]*rs[2];
      fc[0] = (2.0*Dec_Drs - rs[1]*D2ec_Drs2)*Drs/3.0;
    }

  }else{
    static FLOAT fz20[3] = {
      1.709921,                           /* PW */
      1.709920934161365617563962776245,   /* PW (modified) */
      1.709921                            /* OB */
    };

    FLOAT fz, fpz, z4, ec1, alphac;
    FLOAT ectmp, Dec0_Drs, Dec1_Drs, D2ec1_Drs2, Dalphac_Drs, D2alphac_Drs2, Dec_Dz;
    
    fz  =  FZETA(zeta);
    fpz = DFZETA(zeta);
    z4  = POW(zeta, 4);

    /* ec(rs, 1) */
    dp = (fc == NULL) ? NULL : (&D2ec1_Drs2);
    g(func, 1, rs, &ec1, &Dec1_Drs, dp);

    /* -alpha_c(rs) */
    dp = (fc == NULL) ? NULL : (&D2alphac_Drs2);
    g(func, 2, rs, &alphac, &Dalphac_Drs, dp);

    /* what is parametrized is -alpha, so we change signs */
    alphac = -alphac; Dalphac_Drs = -Dalphac_Drs;
    if(dp) D2alphac_Drs2 = -D2alphac_Drs2;
    
    /* save copies that will be needed later */
    Dec0_Drs = Dec_Drs;

    ectmp   =  ec0 + z4*fz*(ec1 - ec0 - alphac/fz20[func]) + fz*alphac/fz20[func];

    Dec_Drs = Dec0_Drs + z4*fz*(Dec1_Drs - Dec0_Drs - Dalphac_Drs/fz20[func]) + fz*Dalphac_Drs/fz20[func];
    Dec_Dz  = (4.0*POW(zeta, 3)*fz + z4*fpz)*(ec1 - ec0 - alphac/fz20[func])
      + fpz*alphac/fz20[func];
    
    if(ec != NULL)
      *ec = ectmp;

    if(vc != NULL){
      vc[0] = ectmp - (rs[1]/3.0)*Dec_Drs - (zeta - 1.0)*Dec_Dz;
      vc[1] = ectmp - (rs[1]/3.0)*Dec_Drs - (zeta + 1.0)*Dec_Dz;
    }    

    if(fc != NULL){
      FLOAT tmp, fppz, D2ec_Dz2, D2ec_DrsDz;
      FLOAT Drs = -(4.0*M_PI/9.0)*rs[2]*rs[2];

      fppz = 0.0;
      if(zeta > -1.0) fppz += POW(1.0 + zeta, -2.0/3.0);
      if(zeta <  1.0) fppz += POW(1.0 - zeta, -2.0/3.0);
      fppz *= 4.0/(9.0*FZETAFACTOR);

      D2ec_Drs2  = D2ec_Drs2 + z4*fz*(D2ec1_Drs2 - D2ec_Drs2 - D2alphac_Drs2/fz20[func])
	+ fz*D2alphac_Drs2/fz20[func];

      D2ec_DrsDz = 4.0*POW(zeta, 3)*fz*(Dec1_Drs - Dec0_Drs - Dalphac_Drs/fz20[func]) + 
	fpz*(z4*(Dec1_Drs - Dec0_Drs) + (1.0 - z4)*Dalphac_Drs/fz20[func]);

      D2ec_Dz2   = 4.0*POW(zeta, 2)*(ec1 - ec0 - alphac/fz20[func])*(3.0*fz + 2.0*zeta*fpz) + 
	fppz*(z4*(ec1 - ec0) + (1.0 - z4)*alphac/fz20[func]);
    
      tmp = (2.0*Dec_Drs - rs[1]*D2ec_Drs2)*Drs/3.0;

      fc[0] = tmp - D2ec_DrsDz*(zeta - 1.0)*Drs 
	+ (zeta - 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta - 1.0)*D2ec_Dz2);
      
      fc[1] = tmp - D2ec_DrsDz*(zeta - 1.0)*Drs 
	+ (zeta + 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta - 1.0)*D2ec_Dz2);

      fc[2] = tmp - D2ec_DrsDz*(zeta + 1.0)*Drs 
	+ (zeta + 1.0)/dens*(D2ec_DrsDz*rs[1]/3.0 + (zeta + 1.0)*D2ec_Dz2);

    }
  }
}


const xc_func_info_type func_info_lda_c_pw = {
  XC_LDA_C_PW,
  XC_CORRELATION,
  "Perdew & Wang",
  XC_FAMILY_LDA,
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,
  NULL,
  lda_c_pw
};

const xc_func_info_type func_info_lda_c_pw_mod = {
  XC_LDA_C_PW_MOD,
  XC_CORRELATION,
  "Perdew & Wang (modified)",
  XC_FAMILY_LDA,
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)\n"
  "Added extra digits to some constants as in the PBE routine\n"
  "http://dft.rutgers.edu/pubs/PBE.asc",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,
  NULL,
  lda_c_pw
};

const xc_func_info_type func_info_lda_c_ob_pw = {
  XC_LDA_C_OB_PW,
  XC_CORRELATION,
  "Ortiz & Ballone (PW parametrization)",
  XC_FAMILY_LDA,
  "G Ortiz and P Ballone, Phys. Rev. B 50, 1391 (1994)\n"
  "G Ortiz and P Ballone, Phys. Rev. B 56, 9970(E) (1997)\n"
  "JP Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,
  NULL,
  lda_c_pw
};
