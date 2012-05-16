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

#define XC_GGA_X_HJS_PBE 525 /* HJS screened exchange PBE version */

typedef struct{
  FLOAT omega;

  const FLOAT *a, *b; /* pointers to the a and b parameters */
} gga_x_hjs_params;

static const FLOAT a_PBE[] = {0.0159941, 0.0852995, -0.160368, 0.152645, -0.0971263, 0.0422061};
static const FLOAT b_PBE[] = {5.33319, -12.4780, 11.0988, -5.11013, 1.71468, -0.610380, 0.307555, -0.0770547, 0.0334840};

static void
gga_x_hjs_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_hjs_params));

  XC(gga_x_hjs_set_params_)(p, 0.2);

  switch(p->info->number){
  case XC_GGA_X_HJS_PBE:
    ((gga_x_hjs_params *)(p->params))->a = a_PBE;
    ((gga_x_hjs_params *)(p->params))->b = b_PBE;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_hjs_init\n");
    exit(1);
  }
}

void 
XC(gga_x_hjs_set_params)(XC(func_type) *p, FLOAT omega)
{
  assert(p != NULL && p->gga != NULL);
  XC(gga_x_hjs_set_params_)(p->gga, omega);
}


void 
XC(gga_x_hjs_set_params_)(XC(gga_type) *p, FLOAT omega)
{
  gga_x_hjs_params *params;

  assert(p->params != NULL);
  params = (gga_x_hjs_params *) (p->params);

  params->omega = omega;
}


#define HEADER 3

/* This implementation follows the one from nwchem */

static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, FLOAT ds,
     FLOAT *f, FLOAT *dfdx, FLOAT *lvrho)
{
  static const FLOAT AA=0.757211, BB=-0.106364, CC=-0.118649, DD=0.609650, EE=-0.0477963;
  static const FLOAT m89=-8.0/9.0;

  FLOAT omega, kF, ss, ss2;
  FLOAT H, F, EG;
  FLOAT nu, zeta, eta, lambda, lambda2, lambda3, lambda4, chi, chi2, chi3, chi4, chi5;
  FLOAT sqzpn2, sqepn2, sqlpn2;
  FLOAT term1, term2, term3, term4, term5, term6;

  FLOAT dnudrho, dssdx, dHds, dFds, dEGds;
  FLOAT dzeta, dchi, dchids, dchidnu;

  assert(p->params != NULL);
  omega = ((gga_x_hjs_params *)(p->params))->omega;

  kF  = POW(3.0*M_PI*M_PI*ds, 1.0/3.0);
  nu  = omega/kF;

  /*  Rescaling the s values to ensure the Lieb-Oxford bound for s>8.3 */
  ss  = X2S*x;
  ss2 = ss*ss; 

  if(order >= 1){
    dnudrho = -nu/(3.0*ds);
    dssdx  = X2S;
  }

  /* first let us calculate H(s) */
  {
    const FLOAT *a, *b;
    FLOAT Hnum, Hden, dHnum, dHden;

    a =  ((gga_x_hjs_params *)(p->params))->a;
    b =  ((gga_x_hjs_params *)(p->params))->b;

    Hnum = ss2*(a[0] + ss*(a[1] + ss*(a[2] + ss*(a[3] + ss*(a[4] + ss*a[5])))));
    Hden = 1.0 + ss*(b[0] + ss*(b[1] + ss*(b[2] + ss*(b[3] + ss*(b[4] + ss*(b[5] + ss*(b[6] + ss*(b[7] + ss*b[8]))))))));

    H = Hnum/Hden;

    if(order >= 1){
      dHnum = ss*(2.0*a[0] + ss*(3.0*a[1] + ss*(4.0*a[2] + ss*(5.0*a[3] + ss*(6.0*a[4] + ss*7.0*a[5])))));
      dHden = b[0] + ss*(2.0*b[1] + ss*(3.0*b[2] + ss*(4.0*b[3] + 
	      ss*(5.0*b[4] + ss*(6.0*b[5] + ss*(7.0*b[6] + ss*(8.0*b[7] + ss*9.0*b[8])))))));

      dHds  = (Hden*dHnum - Hnum*dHden)/(Hden*Hden);
    }
  }

  /* auxiliary variables */
  {
    FLOAT aux, saux;

    zeta   = ss2*H;
    eta    = AA + zeta;
    lambda = DD + zeta;

    aux    = lambda + nu*nu;
    saux   = sqrt(aux);
    chi    = nu/saux;
    
    lambda2 = lambda*lambda;
    lambda3 = lambda*lambda2;
    lambda4 = lambda*lambda3;
    
    chi2 = chi*chi;
    chi3 = chi*chi2;
    chi4 = chi*chi3;
    chi5 = chi*chi4;

    if(order >= 1){
      dzeta   = 2*ss*H + ss2*dHds;
      /* deta = dlambda = dzeta */
      dchids  = -nu*dzeta/(2.0*aux*saux);
      dchidnu = lambda/(aux*saux);
    }
  }

  /* now we calculate F(s) */
  {
    FLOAT aux = 1.0 + 0.25*ss2;

    F = 1.0 - ss2/(27.0*CC*aux) - zeta/(2.0*CC);

    if(order >= 1){
      dFds = -2.0*ss/(27.0*CC*aux*aux) - dzeta/(2.0*CC);
    }
  }
  
  /* and now G(s) */
  {
    FLOAT sqrtl = sqrt(lambda), sqrtz = sqrt(zeta), sqrte = sqrt(eta);

    EG = -(2.0/5.0)*CC*F*lambda - (4.0/15.0)*BB*lambda2 - (6.0/5.0)*AA*lambda3
      - lambda3*sqrtl*((4.0/5.0)*M_SQRTPI + (12.0/5.0)*(sqrtz - sqrte));

    dEGds = -(2.0/5.0)*CC*(dFds*lambda + F*dzeta) - (8.0/15.0)*BB*lambda*dzeta - (18.0/5.0)*AA*lambda2*dzeta
      - (14.0/5.0)*M_SQRTPI*lambda2*sqrtl*dzeta
      - (42.0/5.0)*lambda2*sqrtl*dzeta*((sqrtz - sqrte) + (1.0/7.0)*lambda*(1.0/sqrtz - 1.0/sqrte));
  }

  sqzpn2 = sqrt(zeta + nu*nu);
  sqepn2 = sqrt(eta + nu*nu);
  sqlpn2 = sqrt(lambda + nu*nu);

  term1 = -(4.0/9.0)*BB*(1.0 - chi)/lambda;
  term2 = -(2.0/9.0)*CC*F*(2.0 - 3.0*chi + chi3)/lambda2;
  term3 = -(1.0/9.0)*EG*(8.0 - 15.0*chi + 10.0*chi3 - 3.0*chi5)/lambda3;
  term4 =  2.0*nu*(sqzpn2 - sqepn2);
  term5 =  2.0*zeta*LOG((nu + sqzpn2)/(nu + sqlpn2));
  term6 = -2.0*eta*LOG((nu + sqepn2)/(nu + sqlpn2));

  *f = AA + term1 + term2 + term3 + term4 + term5 + term6;

  if(order >= 1){
    FLOAT dterm1ds, dterm2ds, dterm3ds, dterm4ds, dterm5ds, dterm6ds;
    FLOAT dterm1dnu, dterm2dnu, dterm3dnu, dterm4dnu, dterm5dnu, dterm6dnu;

    dterm1ds = (4.0/9.0)*BB*(lambda*dchids + (1.0 - chi)*dzeta)/lambda2;
    dterm2ds =-(2.0/9.0)*CC*
      (chi - 1.0)*(-dFds*(2.0 - chi - chi2)*lambda + F*(3.0*(1.0 + chi)*lambda*dchids + 2.0*(2.0 - chi - chi2)*dzeta))/lambda3;
    dterm3ds = -(1.0/9.0)*(chi - 1.0)*(chi - 1.0)*
      ((8.0 + chi - 6.0*chi2 - 3.0*chi3)*lambda*dEGds + 
       3.0*EG*(-5.0*(1.0 + chi)*(1.0 + chi)*lambda*dchids + (chi - 1.0)*(8.0 + 9.0*chi + 3.0*chi2)*dzeta))/lambda4;
    dterm4ds = nu*dzeta*(1.0/sqzpn2 - 1.0/sqepn2);
    dterm5ds = dzeta*(-(zeta/lambda)*(1.0 - nu/sqlpn2) + 1.0 + 2.0*LOG((nu + sqzpn2)/(nu + sqlpn2)) - nu/sqzpn2);
    dterm6ds =-dzeta*(-( eta/lambda)*(1.0 - nu/sqlpn2) + 1.0 + 2.0*LOG((nu + sqepn2)/(nu + sqlpn2)) - nu/sqepn2);

    *dfdx = dterm1ds + dterm2ds + dterm3ds + dterm4ds + dterm5ds + dterm6ds;

    dterm1dnu = (4.0/9.0)*BB*dchidnu/lambda;
    dterm2dnu = (2.0/3.0)*CC*F*(1.0 - chi2)*dchidnu/lambda2;
    dterm3dnu = (5.0/3.0)*EG*(1.0 - 2.0*chi2 + chi4)*dchidnu/lambda3;
    dterm4dnu = 2.0*(sqzpn2 - sqepn2 + nu*nu*(1.0/sqzpn2 - 1.0/sqepn2));
    dterm5dnu = 2.0*zeta*(1.0/sqzpn2 - 1.0/sqlpn2);
    dterm6dnu =-2.0*eta*(1.0/sqepn2 - 1.0/sqlpn2);

    *lvrho = dterm1dnu + dterm2dnu + dterm3dnu + dterm4dnu + dterm5dnu + dterm6dnu;
  }

  /* scale and convert to the right variables */
  *dfdx  *= dssdx;
  *lvrho *= dnudrho;
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_hjs_pbe) = {
  XC_GGA_X_HJS_PBE,
  XC_EXCHANGE,
  "HJS screened exchange PBE version",
  XC_FAMILY_GGA,
  "TM Henderson, BG Janesko, and GE Scuseria, J. Chem. Phys. 128, 194105 (2008)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_hjs_init,
  NULL, NULL, 
  work_gga_x
};
