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

#ifndef _LDA_H
#define _LDA_H

#include <math.h>
#include <float.h>
#include "xc_config.h"

/* If strict ANSI, then some useful macros are not defined */
#if defined(__STRICT_ANSI__)
# define M_E            2.7182818284590452354   /* e */
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
double bessi0(double x);
double bessi1(double x);
double bessk0(double x);
double bessk1(double x);
double expint(double x);

/* integration */
typedef void integr_fn(FLOAT *x, int n, void *ex);
FLOAT integrate(integr_fn func, void *ex, FLOAT a, FLOAT b);
void rdqagse(integr_fn f, void *ex, FLOAT *a, FLOAT *b, 
	     FLOAT *epsabs, FLOAT *epsrel, int *limit, FLOAT *result,
	     FLOAT *abserr, int *neval, int *ier, FLOAT *alist__,
	     FLOAT *blist, FLOAT *rlist, FLOAT *elist, int *iord, int *last);
  
#define M_C 137.0359996287515 /* speed of light */

#define RS(x)          (pow((3.0/(4*M_PI*x)), 1.0/3.0))
#define X_FACTOR_C     0.9305257363491000250020102180716672510262     /* 3/8*cur(3/pi)*4^(2/3) */
#define X_FACTOR_2D_C  1.504505556127350098528211870828726895584      /* 8/(3*sqrt(pi))        */
#define X2S            0.1282782438530421943003109254455883701296     /* 1/(2*(6*pi^2)^(1/3))  */
#define FZETAFACTOR    0.519842099789746380
#define FZETA(x)       ((pow(1.0 + (x),  4.0/3.0) + pow(1.0 - (x),  4.0/3.0) - 2.0)/FZETAFACTOR)
#define DFZETA(x)      ((pow(1.0 + (x),  1.0/3.0) - pow(1.0 - (x),  1.0/3.0))*(4.0/3.0)/FZETAFACTOR)
#define D2FZETA(x)     ((4.0/9.0)/FZETAFACTOR)* \
  (ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -2.0/3.0) + pow(1.0 - (x), -2.0/3.0)))
#define D3FZETA(x)     (-(8.0/27.0)/FZETAFACTOR)* \
  (ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -5.0/3.0) - pow(1.0 - (x), -5.0/3.0)))

#define MIN_DENS             1.0e-20
#define MIN_GRAD             1.0e-20
#define MIN_TAU              1.0e-20

#include "xc.h"

void XC(rho2dzeta)(int nspin, const FLOAT *rho, FLOAT *d, FLOAT *zeta);

/* LDAs */
typedef struct XC(lda_rs_zeta) {
  int   order; /* to which order should I return the derivatives */
  FLOAT rs[3], zeta;

  FLOAT zk;
  FLOAT dedrs, dedz;                         /*  first derivatives of zk */
  FLOAT d2edrs2, d2edrsz, d2edz2;            /* second derivatives of zk */
  FLOAT d3edrs3, d3edrs2z, d3edrsz2, d3edz3; /*  third derivatives of zk */
} XC(lda_rs_zeta);

void XC(lda_fxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *fxc);
void XC(lda_kxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *kxc);

/* GGAs */
typedef struct XC(perdew_t) {
  int    nspin;
  FLOAT dens, zeta, gdmt;
  FLOAT ecunif, vcunif[2], fcunif[3];

  FLOAT  rs,  kf,  ks,  phi, t;
  FLOAT drs, dkf, dks, dphi, dt, decunif;

  FLOAT d2rs2, d2rskf, d2rsks, d2rsphi,  d2rst,  d2rsecunif;
  FLOAT        d2kf2,  d2kfks, d2kfphi,  d2kft,  d2kfecunif;
  FLOAT                 d2ks2, d2ksphi,  d2kst,  d2ksecunif;
  FLOAT                         d2phi2, d2phit, d2phiecunif;
  FLOAT                                   d2t2,   d2tecunif;
  FLOAT                                           d2ecunif2;
} XC(perdew_t);

void XC(perdew_params)(const XC(gga_type) *gga_p, const FLOAT *rho, const FLOAT *sigma, int order, XC(perdew_t) *pp);
void XC(perdew_potentials)(XC(perdew_t) *pt, const FLOAT *rho, FLOAT e_gga, int order, 
			   FLOAT *vrho, FLOAT *vsigma, FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
int XC(gga_input_init)(const XC(func_info_type) *info, int nspin, const FLOAT *rho,
		       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
void XC(gga_x_b88_set_params)(XC(gga_type) *p, FLOAT beta);
void XC(gga_c_lyp_set_params)(XC(gga_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d);

void XC(gga_x_pbe_enhance)(const XC(gga_type) *p, int order, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2);


/* hybrid GGAs */
void XC(hyb_gga_alloc)(XC(hyb_gga_type) *p);

/* meta GGAs */
void XC(mgga_x_gvt4_func)(int order, FLOAT x, FLOAT z, FLOAT alpha, const FLOAT *d, 
			  FLOAT *h, FLOAT *dhdx, FLOAT *dhdz);

/* LCAs */
void XC(lca_lch_init)(XC(lca_type) *p);
void XC(lca_omc_init)(XC(lca_type) *p);

void XC(lca_s_lch)(FLOAT rs, FLOAT *s, FLOAT *dsdrs);
void XC(lca_s_omc)(FLOAT rs, FLOAT *s, FLOAT *dsdrs);


#endif
