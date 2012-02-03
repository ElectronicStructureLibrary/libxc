/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
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

#define XC_MGGA_X_TPSS          202 /* Perdew, Tao, Staroverov & Scuseria exchange */

/* some parameters */
static FLOAT b=0.40, c=1.59096, e=1.537, kappa=0.804, mu=0.21951;


/* This is Equation (7) from the paper and its derivatives */
static void 
x_tpss_7(int order, FLOAT p, FLOAT z, 
	 FLOAT *qb, FLOAT *dqbdp, FLOAT *dqbdz, FLOAT *d2qbdp2, FLOAT *d2qbdpz, FLOAT *d2qbdz2)
{
  FLOAT a1, a2, a2_3, a2_5, h1, h2;
  FLOAT alpha, dalphadp, dalphadz, dqbdalpha;
  FLOAT d2alphadpz, d2alphadz2, dqb2dalpha2;

  /* Eq. (8) */
  a1    = (1.0/z - 1.0);
  h1    = 5.0/3.0;

  alpha = h1*a1*p;

  /* Eq. (7) */
  a2    = SQRT(1.0 + b*alpha*(alpha - 1.0));
  h2    = 9.0/20.0;

  *qb   = h2*(alpha - 1.0)/a2 + 2.0*p/3.0;

  if(order < 1) return;              /* And now the derivatives */

  /* Eq. (8) */
  dalphadp = h1*a1;
  dalphadz = -h1*p/(z*z);

  a2_3 = a2*a2*a2;
  dqbdalpha = h2*(1.0 + 0.5*b*(alpha - 1.0))/a2_3;

  *dqbdp = dqbdalpha*dalphadp + 2.0/3.0;
  *dqbdz = dqbdalpha*dalphadz;

  if(order < 2) return;

  /* Eq. (8) */
  d2alphadpz = -h1/(z*z);
  d2alphadz2 = -2.0*dalphadz/z;

  a2_5 = a2*a2*a2_3;
  dqb2dalpha2 = -h2*b/(4.0*a2_5)*(-8.0 + 3.0*b + 12.0*alpha - 7.0*alpha*b +  4*alpha*alpha*b);

  *d2qbdp2 = dqb2dalpha2*dalphadp*dalphadp;
  *d2qbdpz = dqb2dalpha2*dalphadp*dalphadz + dqbdalpha*d2alphadpz;
  *d2qbdz2 = dqb2dalpha2*dalphadz*dalphadz + dqbdalpha*d2alphadz2;
}


/* Equation (10) in all it's glory */
static 
void x_tpss_10(int order, FLOAT p, FLOAT z, 
	       FLOAT *x, FLOAT *dxdp, FLOAT *dxdz, FLOAT *d2xdp2, FLOAT *d2xdpz, FLOAT *d2xdz2)
{
  FLOAT x1, dxdp1, dxdz1, d2xdp1, d2xdpz1, d2xdz1;
  FLOAT aux1, aux2, z2, p2;
  FLOAT qb, dqbdp, dqbdz, d2qbdp2, d2qbdpz, d2qbdz2;
  
  FLOAT a1, a1_2, a2, a3, h3, a4, a5, a6, d1, d1_2;

  /* Equation 7 */
  dqbdp = dqbdz = 0.0;
  x_tpss_7(order, p, z, &qb, &dqbdp, &dqbdz, &d2qbdp2, &d2qbdpz, &d2qbdz2);

  z2   = z*z;
  p2   = p*p;
  aux1 = 10.0/81.0;
  aux2 = 9.0/25.0;
  
  /* first we handle the numerator */
  x1    = 0.0;

  a1 = 1.0 + z2;                               /* first term  */
  a1_2 = a1*a1;
  x1 += (aux1 + c*z2/a1_2)*p;

  a2 = 146.0/2025.0;                           /* second term */
  x1 += a2*qb*qb;

  a3 = SQRT(0.5*(aux2*z2 + p2));               /* third term  */
  h3 = -73.0/405;
  x1 += h3*qb*a3;

  a4 = aux1*aux1/kappa;                        /* forth term  */
  x1 += a4*p2;

  a5 = 2.0*SQRT(e)*aux1*aux2;                  /* fifth term  */
  x1 += a5*z2;

  a6 = e*mu;                                   /* sixth term  */
  x1 += a6*p*p2;

  d1   = 1.0 + SQRT(e)*p;                      /* denominator */
  d1_2 = d1*d1;
  *x   = x1/d1_2;

  if(order < 1) return;              /* the derivatives */

  dxdp1 = dxdz1 = 0.0;

  dxdp1 += aux1 + c*z2/a1_2;                   /* first term  */
  dxdz1 += 2.0*c*z*(1.0 - z2)*p/(a1*a1_2);
  
  dxdp1 += 2.0*a2*qb*dqbdp;                    /* second term */
  dxdz1 += 2.0*a2*qb*dqbdz;
  
  dxdp1 += h3*(a3*dqbdp + 0.5*qb*p/a3);        /* third term  */
  dxdz1 += h3*(a3*dqbdz + 0.5*qb*aux2*z/a3);
  
  dxdp1 += a4*2.0*p;                           /* forth term  */

  dxdz1 += a5*2.0*z;                           /* fifth term  */
  
  dxdp1 += a6*3.0*p2;                          /* sixth term  */
  
  *dxdp = (dxdp1*d1 - 2.0*SQRT(e)*x1)/(d1*d1_2);   /* denominator */
  *dxdz = dxdz1/d1_2;

  if(order < 2) return;

  d2xdp1 = d2xdz1 = d2xdpz1 = 0.0;

  d2xdpz1+= 2.0*c*z*(1.0 - z2)/(a1*a1_2);      /* first term  */
  d2xdz1 += 2.0*c*p*(1.0 - 8.0*z2 + 3.0*z2*z2)/(a1_2*a1_2);

  d2xdp1 += 2.0*a2*(dqbdp*dqbdp + qb*d2qbdp2); /* second term */
  d2xdpz1+= 2.0*a2*(dqbdp*dqbdz + qb*d2qbdpz);
  d2xdz1 += 2.0*a2*(dqbdz*dqbdz + qb*d2qbdz2);

  /*                                              third term */
  d2xdp1 += h3*( aux2*z2*qb  + (p2 + aux2*z2)*(2.0*p*dqbdp            + (p2 + aux2*z2)*d2qbdp2))/(4.0*a3*a3*a3);
  d2xdpz1+= h3*(-aux2*p*z*qb + (p2 + aux2*z2)*(aux2*z*dqbdp + p*dqbdz + (p2 + aux2*z2)*d2qbdpz))/(4.0*a3*a3*a3);
  d2xdz1 += h3*( aux2*p2*qb  + (p2 + aux2*z2)*(2.0*aux2*z*dqbdz       + (p2 + aux2*z2)*d2qbdz2))/(4.0*a3*a3*a3);

  d2xdp1 += a4*2.0;                            /* forth term  */

  d2xdz1 += a5*2.0;                            /* fifth term  */

  d2xdp1 += a6*6.0*p;                          /* sixth term  */

  *d2xdp2 = (6.0*e*x1 + d1*(-4.0*SQRT(e)*dxdp1 + d1*d2xdp1))/(d1_2*d1_2);
  *d2xdpz = (d2xdpz1*d1 - 2.0*SQRT(e)*dxdz1)/(d1*d1_2);
  *d2xdz2 = d2xdz1/d1_2;
}


static void 
func(const XC(mgga_type) *pt, XC(work_mgga_x_params) *r)
{
  FLOAT ss, pp, xx, a1, a1_2, zz;
  FLOAT dxxdp, dxxdz, d2xxdp2, d2xxdpz, d2xxdz2;
  FLOAT dpdx, dzdx, dzdt, d2pdx2, d2zdx2, d2zdxt, d2zdt2;
  
  ss = X2S*r->x;
  pp = ss*ss;

  zz = r->x*r->x/(8.0*r->t);

  /* Eq. 10 */
  x_tpss_10(r->order, pp, zz, &xx, &dxxdp, &dxxdz, &d2xxdp2, &d2xxdpz, &d2xxdz2);

  /* Eq. (5) */
  a1 = kappa/(kappa + xx);
  a1_2 = a1*a1;
  
  r->f = 1.0 + kappa*(1.0 - a1);

  if(r->order < 1) return;

  dpdx = 2.0*ss*X2S;
  dzdx = r->x/(4.0*r->t);
  dzdt = -zz/r->t;

  r->dfdx = a1_2*(dxxdp*dpdx + dxxdz*dzdx);
  r->dfdt = a1_2*dxxdz*dzdt;
  r->dfdu = 0.0;

  if(r->order < 2) return;

  d2pdx2 = 2.0*X2S*X2S;
  d2zdx2 = 1.0/(4.0*r->t);
  d2zdxt = -dzdx/r->t;
  d2zdt2 = -2.0*dzdt/r->t;

  r->d2fdx2 = -2.0*(r->dfdx)*(r->dfdx)/(kappa*a1) +
    a1_2*(d2xxdp2*dpdx*dpdx + 2.0*d2xxdpz*dpdx*dzdx + dxxdp*d2pdx2 + d2xxdz2*dzdx*dzdx + dxxdz*d2zdx2);
  r->d2fdt2 = -2.0*(r->dfdt)*(r->dfdt)/(kappa*a1) +
    a1_2*(d2xxdz2*dzdt*dzdt + dxxdz*d2zdt2);
  r->d2fdxt = -2.0*a1*(r->dfdx)*dxxdz*dzdt/kappa +
    a1_2*(d2xxdpz*dpdx*dzdt + d2xxdz2*dzdx*dzdt + dxxdz*d2zdxt);
}

#include "work_mgga_x.c"

XC(func_info_type) XC(func_info_mgga_x_tpss) = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Tao, Perdew, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  "J Tao, JP Perdew, VN Staroverov, and G Scuseria, Phys. Rev. Lett. 91, 146401 (2003)\n"
  "JP Perdew, J Tao, VN Staroverov, and G Scuseria, J. Chem. Phys. 120, 6898 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
