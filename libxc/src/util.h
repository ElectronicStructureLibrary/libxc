#ifndef _LDA_H
#define _LDA_H

#include <math.h>
#include "xc.h"

#if defined(__STRICT_ANSI__)
# define M_PI           3.14159265358979323846  /* pi */
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

#define   _(is, x)   [3*is + x]
#define  __(i, j)    [2*i + j] 
#define ___(i, j, k) [2*(2*i + j) + k] 

void rho2dzeta(int nspin, double *rho, double *d, double *zeta);

/* LDAs */
void lda_x(void *p, double *rho, double *ex, double *vx, double *fx);

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
