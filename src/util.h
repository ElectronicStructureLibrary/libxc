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
#include "xc.h"

#if defined(__STRICT_ANSI__)
# define M_PI           3.14159265358979323846  /* pi */
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
double asinh(double x);
#endif

#define M_C 137.0359996287515 /* speed of light */
#define min(x,y)  ((x<y) ? x : y)
#define max(x,y)  ((x<y) ? y : x)
#define RS(x)     (pow((3.0/(4*M_PI*x)), 1.0/3.0))
#define X_FACTOR_C  0.9305257363491           /* 3/8*cur(3/pi)*4^(2/3) */
#define FZETAFACTOR 0.519842099789746380
#define FZETA(x)  ((pow(1.0 + zeta, 4.0/3.0) + pow(1.0 - zeta, 4.0/3.0) - 2.0)/FZETAFACTOR)
#define DFZETA(x) ((pow(1.0 + zeta, 1.0/3.0) - pow(1.0 - zeta, 1.0/3.0))*(4.0/3.0)/FZETAFACTOR)

#define MIN_DENS             1.0e-20
#define MIN_GRAD             1.0e-20
#define MIN_TAU              1.0e-20

void rho2dzeta(int nspin, const double *rho, double *d, double *zeta);

/* LDAs */
void xc_lda_fxc_fd(const xc_lda_type *p, const double *rho, double *fxc);
void xc_lda_kxc_fd(const xc_lda_type *p, const double *rho, double *kxc);

/* GGAs */
typedef struct perdew_t {
  int    nspin;
  double dens, zeta, gdmt;
  double ecunif, vcunif[2];

  double  rs,  kf,  ks,  phi, t;
  double drs, dkf, dks, dphi, dt, decunif;
} perdew_t;

void perdew_params(xc_gga_type *gga_p, double *rho, double *sigma, perdew_t *pp);
void perdew_potentials(perdew_t *pt, double *rho, double e_gga, double *vrho, double *vsigma);
void gga_x_b88_set_params(xc_gga_type *p, double beta);
void gga_c_lyp_set_params(xc_gga_type *p, double A, double B, double c, double d);

/* hybrid GGAs */
void xc_hyb_gga_alloc(xc_hyb_gga_type *p);

/* meta-GGAs */
void mgga_x_tpss_init(xc_mgga_type *p);
void mgga_c_tpss_init(xc_mgga_type *p);

void mgga_x_tpss_end(xc_mgga_type *p);
void mgga_c_tpss_end(xc_mgga_type *p);

void mgga_x_tpss(xc_mgga_type *p, double *rho, double *grho, double *tau,
		 double *e, double *dedd, double *dedgd, double *dedtau);
void mgga_c_tpss(xc_mgga_type *p, double *rho, double *grho, double *tau,
		 double *e, double *dedd, double *dedgd, double *dedtau);

/* LCAs */
void lca_lch_init(xc_lca_type *p);
void lca_omc_init(xc_lca_type *p);

void lca_s_lch(double rs, double *s, double *dsdrs);
void lca_s_omc(double rs, double *s, double *dsdrs);


#endif
