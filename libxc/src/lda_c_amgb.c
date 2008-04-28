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
 Correlation energy per particle and potentials for a homogeneous electron
 gas in 2D, as parametrized by Attacalite et al.
************************************************************************/

#define XC_LDA_C_AMGB  15   /* Attacalite et al             */


/* parameters necessary to the calculation */
static FLOAT a[3] = { -0.1925,     0.117331,    0.0234188 };
static FLOAT b[3] = {  0.0863136, -3.394e-2,   -0.037093  };
static FLOAT c[3] = {  0.0572384, -7.66765e-3,  0.0163618 };
static FLOAT d[3] = {  0.0,        0.0,         0.0       };
static FLOAT e[3] = {  1.0022,     0.4133,      1.424301  };
static FLOAT f[3] = { -0.02069,    0.0,         0.0       };
static FLOAT g[3] = {  0.33997,    6.68467e-2,  0.0       };
static FLOAT h[3] = {  1.747e-2,   7.799e-4,    1.163099  };
static FLOAT beta = 1.3386, ax = 0.0;


/* Initialization */
static void
lda_c_amgb_init(void *p)
{
  int i;
  
  /* initialize a couple of constants */
  for(i=0; i<3; i++) d[i] = -a[i]*h[i];
  ax = -4.0/(3.0*M_PI*sqrt(2.0));
}


/* Equation [1].4 */
static FLOAT
alpha(int i, FLOAT *rs)
{
  return a[i] + (b[i]*rs[1] + c[i]*rs[2] + d[i]*rs[3]) *
    log(1.0 + 1.0/(e[i]*rs[1] + f[i]*rs[0]*rs[1] + g[i]*rs[2] + h[i]*rs[3]));
}


/* Equation [1].C3 */
static FLOAT 
dalphadrs(int i, FLOAT *rs)
{
  FLOAT efe, efep, lg, x;
  
  efe  = e[i]*rs[1] + f[i]*rs[0]*rs[1] + g[i]*rs[2] + h[i]*rs[3]; /* Eq. [2] C5 */
  efep = e[i] + 1.5*f[i]*rs[0] + 2.0*g[i]*rs[1] + 3.0*h[i]*rs[2]; /* Eq. [2] C6 */
  lg = log(1.0 + 1.0/efe);
  x  = ((b[i]*rs[1] + c[i]*rs[2] + d[i]*rs[3])*efep)/(efe*efe + efe);
  return (b[i] + 2.0*c[i]*rs[1] + 3.0*d[i]*rs[2])*lg - x;         /* Eq. [2] C3 */
}


static void
lda_c_amgb(const void *p_, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *fc)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;
  
  FLOAT dens, zeta, rs[4], eps_c;
  
  assert(p  != NULL);
  assert(ax != 0.0);
  
  /* get the trace and the polarization of the density */
  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);
  
  /* Wigner radius: formula for 2D, of course */
  rs[1] = sqrt(1.0/(M_PI*dens));
  rs[0] = sqrt(rs[1]);
  rs[2] = rs[1]*rs[1];
  rs[3] = rs[2]*rs[1];
  
  if(p->nspin == XC_UNPOLARIZED){
    /* In unpolarized cases, expressions are fairly simple. */
    
    eps_c = alpha(0, rs);

    if(zk != NULL)
      *zk = eps_c;

    if(vrho != NULL)
      vrho[0] = eps_c - 0.5*rs[1]*dalphadrs(0, rs);
    
  }else{ /* XC_POLARIZED */
    /* In the case of spin-polarized calculations, all this is necessary... */
    
    FLOAT ex, ex0, ex6, calf, calfp, decdrs, decdz, zeta2, zeta4;
    
    /* Unpolarized exchange energy */
    ex0 = ((-4.0*sqrt(2.0))/(3.0*M_PI*rs[1]));
    
    /* Polarized exchange energy */
    ex  = 0.5*(POW(1.0 + zeta, 1.5) + POW(1.0 - zeta, 1.5))*ex0;
    
    /* Taylor expansion of ex, in zeta, beyond fourth order. */
    zeta2 = zeta*zeta;
    zeta4 = zeta2*zeta2;
    ex6 = ex - (1.0 + (3.0/8.0)*zeta2 + (3.0/128.0)*zeta4)*ex0;

    /* Correlation energy, Eq. [1] 3 */
    eps_c = alpha(0, rs) + alpha(1, rs)*zeta2 + alpha(2, rs)*zeta4 + (exp(-beta*rs[1]) - 1.0)*ex6;
    
    if(zk != NULL)
      *zk = eps_c;

    if(vrho == NULL) return;

    /* Function calf, Eq. [2] 4.10 */
    calf = POW(1.0 + zeta, 1.5) + POW(1.0 - zeta, 1.5) - 
      (2.0 + (3.0/4.0)*zeta2 + (3.0/64.0)*zeta4);
    
    /* Function calfp, Eq. [2] C8 */
    calfp = 1.5*(sqrt(1.0 + zeta) - sqrt(1.0 - zeta)) 
      - 1.5*zeta - (3.0/16.0)*zeta*zeta2;
    
    /* Derivative of the correlation energy with respect to rs, Eq. [2] C2 */
    decdrs = ax*calf*(1.0 - exp(-beta*rs[1])*(1.0 + beta*rs[1]))/rs[2] + 
      dalphadrs(0, rs) + dalphadrs(1, rs)*zeta2 + dalphadrs(2, rs)*zeta4;
    
    /* Derivative of the correlation energy with respect to zeta, Eq. [2] C7 */
    decdz = ax*(exp(-beta*rs[1]) - 1.0)*calfp/rs[1] + 
      2.0*alpha(1, rs)*zeta + 4.0*alpha(2, rs)*zeta*zeta2;
    
    /* And finally, the potentials, Eq. [2] C1 */
    if(vrho != NULL){
      vrho[0] = eps_c - 0.5*rs[1]*decdrs - (zeta - 1.0)*decdz;
      vrho[1] = eps_c - 0.5*rs[1]*decdrs - (zeta + 1.0)*decdz;
    }
  }
}


const XC(func_info_type) XC(func_info_lda_c_amgb) = {
  XC_LDA_C_AMGB,
  XC_CORRELATION,
  "AMGB (for 2D systems)",
  XC_FAMILY_LDA,
  "C Attacalite et al, Phys. Rev. Lett. 88, 256601 (2002)\n"
  "C Attacalite, PhD thesis",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_c_amgb_init,
  NULL,
  lda_c_amgb
};
