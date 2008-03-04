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
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

#define XC_MGGA_X_TPSS          201 /* Perdew, Tao, Staroverov & Scuseria exchange */

static XC(func_info_type) func_info_mgga_x_tpss = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Perdew, Tao, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  "J.P.Perdew, Tao, Staroverov, and Scuseria, Phys. Rev. Lett. 91, 146401 (2003)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC
};


void XC(mgga_x_tpss_init)(XC(mgga_type) *p)
{
  p->info = &func_info_mgga_x_tpss;
  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_x_init)(p->lda_aux, XC_UNPOLARIZED, 3, XC_NON_RELATIVISTIC);
}


void XC(mgga_x_tpss_end)(XC(mgga_type) *p)
{
  free(p->lda_aux);
}


/* some parameters */
static FLOAT b=0.40, c=1.59096, e=1.537, kappa=0.804, mu=0.21951;


/* This is Equation (7) from the paper and its derivatives */
static void 
x_tpss_7(FLOAT p, FLOAT z, 
	 FLOAT *qb, FLOAT *dqbdp, FLOAT *dqbdz)
{
  FLOAT alpha, dalphadp, dalphadz;

  { /* Eq. (8) */
    FLOAT a = (1.0/z - 1.0), h = 5.0/3.0;
    alpha    = h*a*p;
    dalphadp = h*a;
    dalphadz = -h*p/(z*z);
  }

  { /* Eq. (7) */
    FLOAT dqbda;
    FLOAT a = sqrt(1.0 + b*alpha*(alpha-1.0)), h = 9.0/20.0;
    dqbda = h*(1.0 + 0.5*b*(alpha-1.0))/POW(a, 3);

    *qb    = h*(alpha - 1.0)/a + 2.0*p/3.0;
    *dqbdp = dqbda*dalphadp + 2.0/3.0;
    *dqbdz = dqbda*dalphadz;
  }

}

/* Equation (10) in all it's glory */
static 
void x_tpss_10(FLOAT p, FLOAT z, 
	       FLOAT *x, FLOAT *dxdp, FLOAT *dxdz)
{
  FLOAT x1, dxdp1, dxdz1;
  FLOAT aux1, z2, p2;
  FLOAT qb, dqbdp, dqbdz;
  
  /* Equation 7 */
  x_tpss_7(p, z, &qb, &dqbdp, &dqbdz);

  z2   = z*z;
  p2   = p*p; 
  aux1 = 10.0/81.0;
  
  /* first we handle the numerator */
  x1    = 0.0;
  dxdp1 = 0.0;
  dxdz1 = 0.0;

  { /* first term */
    FLOAT a = 1.0+z2, a2 = a*a;
    x1    += (aux1 + c*z2/a2)*p;
    dxdp1 += (aux1 + c*z2/a2);
    dxdz1 += c*2.0*z*(1.0-z2)*p/(a*a2);
  }
  
  { /* second term */
    FLOAT a = 146.0/2025.0*qb;
    x1    += a*qb;
    dxdp1 += 2.0*a*dqbdp;
    dxdz1 += 2.0*a*dqbdz;
  }
  
  { /* third term */
    FLOAT a = sqrt(0.5*(9.0*z2/25.0 + p2));
    FLOAT h = 72.0/405;
    x1    += -h*qb*a;
    dxdp1 += -h*(a*dqbdp + 0.5*qb*p/a);
    dxdz1 += -h*(a*dqbdz + 0.5*qb*(9.0/25.0)*z/a);
  }
  
  { /* forth term */
    FLOAT a = aux1*aux1/kappa;
    x1    += a*p2;
    dxdp1 += a*2.0*p;
  }
  
  { /* fifth term */
    FLOAT a = 2.0*sqrt(e)*aux1*9.0/25.0;
    x1    += a*z2;
    dxdz1 += a*2.0*z;
  }
  
  { /* sixth term */
    FLOAT a = e*mu;
    x1    += a*p*p2;
    dxdp1 += a*3.0*p2;
  }
  
  /* and now the denominator */
  {
    FLOAT a = 1.0+sqrt(e)*p, a2 = a*a;
    *x    = x1/a2;
    *dxdp = (dxdp1*a - 2.0*sqrt(e)*x1)/(a2*a);
    *dxdz = dxdz1/a2;
  }
}

static void 
x_tpss_para(XC(mgga_type) *pt, FLOAT rho, FLOAT *grho, FLOAT tau_,
	    FLOAT *energy, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau)
{

  FLOAT gdms, p, tau, tauw, z;
  FLOAT x, dxdp, dxdz, Fx, dFxdx;
  FLOAT exunif, vxunif;
  
  tau = max(tau_, MIN_TAU);

  /* get the uniform gas energy and potential */
  XC(lda_vxc)(pt->lda_aux, &rho, &exunif, &vxunif);

  /* calculate |nabla rho|^2 */
  gdms = grho[0]*grho[0] + grho[1]*grho[1] + grho[2]*grho[2];
  gdms = max(MIN_GRAD*MIN_GRAD, gdms);
  
  /* Eq. (4) */
  p = gdms/(4.0*POW(3*M_PI*M_PI, 2.0/3.0)*POW(rho, 8.0/3.0));

  /* von Weisaecker kinetic energy density */
  tauw = gdms/(8.0*rho);
  z  = tauw/tau;

  /* Eq. 10 */
  x_tpss_10(p, z, &x, &dxdp, &dxdz);

  { /* Eq. (5) */
    FLOAT a = kappa/(kappa + x);
    Fx    = 1.0 + kappa*(1.0 - a);
    dFxdx = a*a;
  }
  
  { /* Eq. (3) */
    int i;
    FLOAT a = rho*exunif*dFxdx;

    *energy = exunif*Fx;
    *dedd   = vxunif*Fx + exunif*dFxdx*(-(8.0/3.0)*p*dxdp - z*dxdz);
    *dedtau = a * (-z/tau*dxdz);

    for(i=0; i<3; i++)
      dedgd[i] = a * 2.0*grho[i]/gdms * (p*dxdp + z*dxdz);
  }
}


void 
XC(mgga_x_tpss)(XC(mgga_type) *p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
	    FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau)
{
  if(p->nspin == XC_UNPOLARIZED){
    x_tpss_para(p, rho[0], grho, tau[0], e, dedd, dedgd, dedtau);

  }else{ 
    /* The spin polarized version is handle using the exact spin scaling
          Ex[n1, n2] = (Ex[2*n1] + Ex[2*n2])/2
    */
    int is;

    *e = 0.0;
    for(is=0; is<2; is++){
      /* FLOAT gr[3], e1;
         int i;
         for(i=0; i<3; i++) gr[i] = 2.0*grho _(is, i);

         x_tpss_para(p, 2.0*rho[is], gr, 2.0*tau[is], &e1, 
      	 &(dedd[is]), &(dedgd _(is, 0)), &(dedtau[is]));
	 *e += e1; */
    }
  }
}
