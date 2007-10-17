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
#include <math.h>
#include "util.h"

#define XC_GGA_C_PW91 134 /* Perdew & Wang 91 */

static double pw91_nu, pw91_beta;
static const double
  pw91_C_c0  = 4.235e-3, 
  pw91_alpha = 0.09;

static void gga_c_pw91_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_init(p->lda_aux, XC_LDA_C_PW_MOD, p->nspin);

  pw91_nu   = 16.0/M_PI * pow(3.0*M_PI*M_PI, 1.0/3.0);
  pw91_beta = pw91_nu*pw91_C_c0;
}


static void gga_c_pw91_end(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  free(p->lda_aux);
}


inline void
A_eq14(double ec, double g, double *A, double *dec, double *dg)
{
  double g2, g3, dd;

  g2 = g*g;
  g3 = g*g2;

  dd = -2.0*pw91_alpha*ec/(g3*pw91_beta*pw91_beta);
  dd = exp(dd);

  *A   = (2.0*pw91_alpha/pw91_beta) / (dd - 1.0);

  *dec = -(*A)/(dd - 1.0) * dd * (-2.0*pw91_alpha/(g3*pw91_beta*pw91_beta));
  *dg  = -(*A)/(dd - 1.0) * dd * (+6.0*pw91_alpha*ec/(g*g3*pw91_beta*pw91_beta));
}

inline void
H0_eq13(double   ec, double   g, double   t, double *H0,
	double *dec, double *dg, double *dt)
{
  double A, dAdec, dAdg;
  double g3, t2, t4, n0, d0, dd, dA;

  A_eq14(ec, g, &A, &dAdec, &dAdg);

  g3 = g*g*g;
  t2 = t*t;
  t4 = t2*t2;

  n0 = t2 + A*t4;
  d0 = 1.0 + A*t2 + A*A*t4;

  *H0 = g3 * pw91_beta*pw91_beta/(2.0*pw91_alpha) *
    log(1.0 + 2.0*pw91_alpha/pw91_beta * n0/d0);

  dd = d0*(pw91_beta + n0*(2.0*pw91_alpha + A*pw91_beta));
  dA = -A*g3*t2*t4*(2.0 + A*t2)*pw91_beta*pw91_beta/dd;

  *dec = dA * dAdec;
  *dg  = (*H0)*3.0/g + dA*dAdg;
  *dt  = 2.0*g3*t*(1.0 + 2.0*A*t2)*pw91_beta*pw91_beta/dd;
}


/* pade parametrized form of C-xc found in
   M Rasolt & DJW Geldart, Phys. Rev. B 34, 1325 (1986)
*/
inline void 
Rasold_Geldart_C_xc(double rs, double *C_xc, double *drs)
{
  const double 
    a[3] = {2.568, 23.266, 0.007389},
    b[3] = {1.0, 8.723, 0.472};
  
  double d0, d1, n0, n1;

  n0 = (a[0] + rs*(a[1] + rs*a[2]));
  d0 =  b[0] + rs*(b[1] + rs*(b[2] + 10.0*rs*a[2]));

  n1 = a[1] + 2.0*rs*a[2];
  d1 = b[1] + 2.0*rs*b[2] + 10.0*3.0*rs*rs*a[2];

  *C_xc = n0/(1000.0*d0);
  *drs  = (n1*d0 - n0*d1)/(1000.0*d0*d0);
}


inline void 
H1_eq15(double   rs, double   g, double   t, double   ks, double   kf, double *H1,
	double *drs, double *dg, double *dt, double *dks, double *dkf)
{
  const double C_xc0 = 2.568e-3, C_x = -0.001667212;

  double g3, g4, t2, kf2, ks2, dd1, dd2;
  double C_xc, dC_xc;

  g3  = g*g*g;
  g4  = g3*g;
  t2  = t*t;
  ks2 = ks*ks;
  kf2 = kf*kf;

  dd1 = -100.0 * g4 * (ks2/kf2) * t2;
  dd1 = exp(dd1);

  Rasold_Geldart_C_xc(rs, &C_xc, &dC_xc);
  dd2 = C_xc - C_xc0 - 3.0*C_x/7.0;

  *H1  = pw91_nu * dd2 * g3 * t2 * dd1;

  *drs = pw91_nu * dC_xc * g3 * t2 * dd1;
  *dg  = (*H1) * (3.0/g - 100.0* 4.0*g3 *(ks2/kf2)*t2); /* g can not be zero */
  *dt  = pw91_nu * dd2 * g3 * dd1 * 
    (2.0*t - t2*100.0*g4*(ks2/kf2)* 2.0*t);
  *dks = (*H1) * (-100.0*g4*( 2.0*ks /kf2)*t2);
  *dkf = (*H1) * (+100.0*g4*( 2.0*ks2/(kf*kf2) )*t2);
}


inline void 
ec_eq9(double   ec, double   rs, double   t, double   g, double   ks, double   kf, double  *ec_gga,
       double *dec, double *drs, double *dt, double *dg, double *dks, double *dkf)
{
  double H0, dH0dec, dH0dg, dH0dt;
  double H1, dH1drs, dH1dg, dH1dt, dH1dks, dH1dkf;

  H0_eq13(ec, g, t, &H0, 
	  &dH0dec, &dH0dg, &dH0dt);
  H1_eq15(rs, g, t, ks, kf, &H1,
	  &dH1drs, &dH1dg, &dH1dt, &dH1dks, &dH1dkf);

  *ec_gga = ec + H0 + H1;
  *dec    = 1.0 + dH0dec;
  *drs    = dH1drs;
  *dt     = dH0dt + dH1dt;
  *dg     = dH0dg + dH1dg;
  *dks    = dH1dks;
  *dkf    = dH1dkf;
}

static void gga_c_pw91(void *p_, double *rho, double *sigma,
		       double *e, double *vrho, double *vsigma)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  double dens, zeta, ecunif, vcunif[2];
  double rs, kf, ks, phi, gdmt, t;

  xc_lda_vxc(p->lda_aux, rho, &ecunif, vcunif);
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs = RS(dens);
  kf = pow(3.0*M_PI*M_PI*dens, 1.0/3.0);
  ks = sqrt(4.0*kf/M_PI);

  phi  = 0.5*(pow(1.0 + zeta, 2.0/3.0) + pow(1.0 - zeta, 2.0/3.0));

  /* get gdmt = |nabla n| */
  gdmt = sigma[0];
  if(p->nspin == XC_POLARIZED) gdmt += 2.0*sigma[1] + sigma[2];
  gdmt = sqrt(gdmt);
  if(gdmt < MIN_GRAD) gdmt = MIN_GRAD;

  t = gdmt/(2.0*phi*ks*dens);

  {
    double e_gga, dec, drs, dt, dphi, dks, dkf;
    double drsdd, dkfdd, dksdd, dzdd[2], dpdz;
    int is;

    ec_eq9(ecunif, rs, t, phi, ks, kf, &e_gga,
	   &dec, &drs, &dt, &dphi, &dks, &dkf);
    //fprintf(stderr, "%lf %.10lf %.10lf\n", ecunif, e_gga, dec);

    *e = e_gga;

    drsdd   = -rs/(3.0*dens);
    dkfdd   =  kf/(3.0*dens);
    dksdd   = 0.5*ks*dkfdd/kf;
    dzdd[0] =  (1.0 - zeta)/dens;
    dzdd[1] = -(1.0 + zeta)/dens;
    dpdz    = 0.0;
    if(fabs(1.0 + zeta) >= MIN_DENS)
      dpdz += (1.0/3.0)/pow(1.0 + zeta, 1.0/3.0);
    if(fabs(1.0 - zeta) >= MIN_DENS)
      dpdz -= (1.0/3.0)/pow(1.0 - zeta, 1.0/3.0);
    
    /* add the t contributions to the other derivatives */
    dphi += dt * (-t/phi);
    dks  += dt * (-t/ks);

    /* calculate vrho */
    for(is=0; is<p->nspin; is++){
      if(rho[is] > MIN_DENS){
	double decudd;

	vrho[is]  = e_gga;

	decudd = (vcunif[is] - ecunif)/dens;
	vrho[is] += dens * (dec*decudd + drs*drsdd + dkf*dkfdd + dks*dksdd - dt*t/dens);
	vrho[is] += dens * dphi*dpdz*dzdd[is];
      }else{
	vrho[is] = 0.0;
      }
    }
    
    { /* calculate now vsigma */
      double dtdsig;

      dtdsig  = t/(2.0*gdmt*gdmt);
      vsigma[0] = dens*dt*dtdsig;
      if(is == 2){
	vsigma[1] = 2.0*vsigma[0];
	vsigma[2] =     vsigma[0];
      }
    }
  }
}

const xc_func_info_type func_info_gga_c_pw91 = {
  XC_GGA_C_PW91,
  XC_EXCHANGE,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 46, 6671 (1992)\n"
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 48, 4978(E) (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_c_pw91_init,
  gga_c_pw91_end,
  NULL,            /* this is not an LDA                   */
  gga_c_pw91,
};
