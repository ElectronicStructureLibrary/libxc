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
#include <assert.h>
#include "util.h"


/************************************************************************
   L. Hedin and  B.I. Lundqvist
   O. Gunnarsson and B. I. Lundqvist
************************************************************************/

#define XC_LDA_C_HL  4   /* Hedin & Lundqvist            */
#define XC_LDA_C_GL  5   /* Gunnarson & Lundqvist        */

static void hl_f(int func, int i, FLOAT rs, FLOAT *ec, FLOAT *vc)
{
  static const 
    FLOAT r[2][2] = {{21.0,   21.0},     /* HL unpolarized only*/
		      {11.4,   15.9}};    /* GL */
  static const 
    FLOAT c[2][2] = {{ 0.0225, 0.0225},  /* HL unpolarized only */
		      { 0.0333, 0.0203}}; /* GL */
  
  FLOAT a, x, x2, x3;
  
  x   = rs/r[func][i];
  x2  = x*x;
  x3  = x2*x;
  
  a   = log(1.0 + 1.0/x);
  *ec = -c[func][i]*((1.0 + x3)*a - x2 + 0.5*x - 1.0/3.0);
  
  *vc = -c[func][i]*a;
}


static void lda_c_hl(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  FLOAT ecp, vcp;
  FLOAT dens, zeta, rs;
  int func = p->info->number - XC_LDA_C_HL;

  /* sanity check */
  assert(func==0 || func==1);
  
  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);
  rs = RS(dens); /* Wigner radius */  

  hl_f(func, 0, rs, &ecp, &vcp);
  
  if(p->nspin==XC_UNPOLARIZED){
    *ec   = ecp;
    vc[0] = vcp;
    
  }else{ /* XC_POLARIZED */
    FLOAT ecf, vcf, fz, dfz, t1;
    
    fz  =  FZETA(zeta);
    dfz = DFZETA(zeta);
    
    hl_f(func, 1, rs, &ecf, &vcf);
    
    *ec = ecp + (ecf - ecp)*fz;                  /* the energy    */
    t1  = vcp + (vcf - vcp)*fz;
    
    vc[0] = t1 + (ecf - ecp)*dfz*( 1.0 - zeta);  /* the potential */
    vc[1] = t1 + (ecf - ecp)*dfz*(-1.0 - zeta);
  }
}

const XC(func_info_type) XC(func_info_lda_c_hl) = {
  XC_LDA_C_HL,
  XC_CORRELATION,
  "Hedin & Lundqvist",
  XC_FAMILY_LDA,
  /* can someone get me this paper, so I can find all coefficients? */
  "L Hedin and BI Lundqvist, J. Phys. C 4, 2064 (1971)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_hl,     /* lda  */
};

const XC(func_info_type) XC(func_info_lda_c_gl) = {
  XC_LDA_C_GL,
  XC_CORRELATION,
  "Gunnarson & Lundqvist",
  XC_FAMILY_LDA,
  "O Gunnarsson and BI Lundqvist, PRB 13, 4274 (1976)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_hl,     /* lda  */
};
