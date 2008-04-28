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

static void hl_f(int func, int i, FLOAT rs, FLOAT *zk, FLOAT *drs, FLOAT *d2rs)
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
  *zk = -c[func][i]*((1.0 + x3)*a - x2 + 0.5*x - 1.0/3.0);
  
  if(drs != NULL)
    *drs = -c[func][i]/r[func][i]*(3.0*x*(x*a - 1) - 1/x + 3.0/2.0);

  if(d2rs != NULL)
    *d2rs = -c[func][i]/(r[func][i]*r[func][i]*x2*(1.0 + x))*
      (1.0 + x - 3.0*x2 - 6.0*x3*(1.0 - (1.0 + x)*a));
}


static inline void 
func(const XC(lda_type) *p, FLOAT *rs, FLOAT zeta, 
     FLOAT *zk, FLOAT *dedrs, FLOAT *dedz, 
     FLOAT *d2edrs2, FLOAT *d2edrsz, FLOAT *d2edz2)
{
  int func = p->info->number - XC_LDA_C_HL;
  assert(func==0 || func==1);
  
  hl_f(func, 0, rs[1], zk, dedrs, d2edrs2);

  if(p->nspin==XC_UNPOLARIZED) return; /* nothing else to do */

  { /* XC_POLARIZED */
    FLOAT ecp, vcp, fcp;
    FLOAT ecf, vcf, fcf, fz, dfz, d2fz;
    
    /* store paramagnetic values */
    ecp = *zk;
    if(dedrs   != NULL) vcp = *dedrs;
    if(d2edrs2 != NULL) fcp = *d2edrs2;

    /* get ferromagnetic values */
    hl_f(func, 1, rs[1], &ecf, dedrs, d2edrs2);
    if(dedrs   != NULL) vcf = *dedrs;
    if(d2edrs2 != NULL) fcf = *d2edrs2;
    
    fz  =  FZETA(zeta);
    *zk = ecp + (ecf - ecp)*fz;

    if(dedrs==NULL && d2edrs2==NULL) return; /* nothing else to do */

    dfz = DFZETA(zeta);
    if(dedrs!=NULL){
      *dedrs = vcp + (vcf - vcp)*fz;
      *dedz  = (ecf - ecp)*dfz;
    }

    if(d2edrs2==NULL) return; /* nothing else to do */
    
    d2fz = D2FZETA(zeta);
    *d2edrs2 = fcp + (fcf - fcp)*fz;
    *d2edrsz =       (vcf - vcp)*dfz;
    *d2edz2  =       (ecf - ecp)*d2fz;
  }
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_hl) = {
  XC_LDA_C_HL,
  XC_CORRELATION,
  "Hedin & Lundqvist",
  XC_FAMILY_LDA,
  /* can someone get me this paper, so I can find all coefficients? */
  "L Hedin and BI Lundqvist, J. Phys. C 4, 2064 (1971)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};

const XC(func_info_type) XC(func_info_lda_c_gl) = {
  XC_LDA_C_GL,
  XC_CORRELATION,
  "Gunnarson & Lundqvist",
  XC_FAMILY_LDA,
  "O Gunnarsson and BI Lundqvist, PRB 13, 4274 (1976)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
