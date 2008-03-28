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

#ifndef _LDA_H
#define _LDA_H

#include <math.h>
#include "xc_config.h"

/* If strict ANSI, then some useful macros are not defined */
#if defined(__STRICT_ANSI__)
# define M_PI           3.14159265358979323846  /* pi */
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
double asinh (double x);
float  asinhf(float  x);
#endif

/* Very useful macros */
#define min(x,y)  ((x<y) ? x : y)
#define max(x,y)  ((x<y) ? y : x)

/* special functions */
double lambert_w(double z);

#define M_C 137.0359996287515 /* speed of light */

#define RS(x)     (pow((3.0/(4*M_PI*x)), 1.0/3.0))
#define X_FACTOR_C  0.9305257363491           /* 3/8*cur(3/pi)*4^(2/3) */
#define FZETAFACTOR 0.519842099789746380
#define FZETA(x)  ((pow(1.0 + zeta, 4.0/3.0) + pow(1.0 - zeta, 4.0/3.0) - 2.0)/FZETAFACTOR)
#define DFZETA(x) ((pow(1.0 + zeta, 1.0/3.0) - pow(1.0 - zeta, 1.0/3.0))*(4.0/3.0)/FZETAFACTOR)

#define MIN_DENS             1.0e-20
#define MIN_GRAD             1.0e-20
#define MIN_TAU              1.0e-20

#include "xc.h"

void XC(rho2dzeta)(int nspin, const FLOAT *rho, FLOAT *d, FLOAT *zeta);

/* LDAs */
void XC(lda_fxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *fxc);
void XC(lda_kxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *kxc);

/* GGAs */
typedef struct XC(perdew_t) {
  int    nspin;
  FLOAT dens, zeta, gdmt;
  FLOAT ecunif, vcunif[2];

  FLOAT  rs,  kf,  ks,  phi, t;
  FLOAT drs, dkf, dks, dphi, dt, decunif;
} XC(perdew_t);

void XC(perdew_params)(XC(gga_type) *gga_p, FLOAT *rho, FLOAT *sigma, XC(perdew_t) *pp);
void XC(perdew_potentials)(XC(perdew_t) *pt, FLOAT *rho, FLOAT e_gga, FLOAT *vrho, FLOAT *vsigma);
void XC(gga_x_b88_set_params)(XC(gga_type) *p, FLOAT beta);
void XC(gga_c_lyp_set_params)(XC(gga_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d);

/* hybrid GGAs */
void XC(hyb_gga_alloc)(XC(hyb_gga_type) *p);

/* meta-GGAs */
void XC(mgga_x_tpss_init)(XC(mgga_type) *p);
void XC(mgga_c_tpss_init)(XC(mgga_type) *p);

void XC(mgga_x_tpss_end)(XC(mgga_type) *p);
void XC(mgga_c_tpss_end)(XC(mgga_type) *p);

void XC(mgga_x_tpss)(XC(mgga_type) *p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
		 FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau);
void XC(mgga_c_tpss)(XC(mgga_type) *p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
		 FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau);

/* LCAs */
void XC(lca_lch_init)(XC(lca_type) *p);
void XC(lca_omc_init)(XC(lca_type) *p);

void XC(lca_s_lch)(FLOAT rs, FLOAT *s, FLOAT *dsdrs);
void XC(lca_s_omc)(FLOAT rs, FLOAT *s, FLOAT *dsdrs);


#endif
