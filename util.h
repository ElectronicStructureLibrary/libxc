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
#define FZETAFACTOR 0.519842099789746380
#define FZETA(x)  ((pow(1.0 + zeta, 4.0/3.0) + pow(1.0 - zeta, 4.0/3.0) - 2.0)/FZETAFACTOR)
#define DFZETA(x) ((pow(1.0 + zeta, 1.0/3.0) - pow(1.0 - zeta, 1.0/3.0))*(4.0/3.0)/FZETAFACTOR)

#define MIN_DENS             1.0e-14
#define MIN_GRAD             1.0e-14
#define MIN_TAU              1.0e-14

#define   _(is, x)   [3*is + x]
#define  __(i, j)    [2*i + j] 
#define ___(i, j, k) [2*(2*i + j) + k] 

void rho2dzeta(int nspin, double *rho, double *d, double *zeta);

/* LDAs */
void lda_c_wigner_init (lda_type *p);
void lda_c_rpa_init    (lda_type *p);
void lda_c_hl_init     (lda_type *p);
void lda_c_gl_init     (lda_type *p);
void lda_c_vwn_init    (lda_type *p);
void lda_c_vwn_rpa_init(lda_type *p);
void lda_c_pz_init     (lda_type *p);
void lda_c_pz_mod_init (lda_type *p);
void lda_c_ob_pz_init  (lda_type *p);
void lda_c_pw_init     (lda_type *p);
void lda_c_ob_pw_init  (lda_type *p);
void lda_c_lyp_init    (lda_type *p);
void lda_c_amgb_init   (lda_type *p);

void lda_x       (lda_type *p, double *rho, double *ex, double *vx, double *fx);
void lda_c_wigner(lda_type *p, double rs, double *ec, double *vc);
void lda_c_rpa   (lda_type *p, double rs, double *ec, double *vc);
void lda_c_hl    (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_xalpha(lda_type *p, double *rho, double *ec, double *vc, double *fx);
void lda_c_vwn   (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_pz    (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_pw    (lda_type *p, double rs, double dens, double zeta, double *ec, double *vc, double *fx);
void lda_c_lyp   (lda_type *p, double rs, double *ec, double *vc);
void lda_c_amgb  (lda_type *p, double *rho, double *ec, double *vc);

/* GGAs */
void gga_x_pbe_init(gga_type *p);
void gga_c_pbe_init(gga_type *p);

void gga_pbe_end(gga_type *p);
void gga_lb_end (gga_type *p);

void gga_x_pbe(gga_type *p, double *rho, double *grho,
	       double *e, double *dedd, double *dedgd);
void gga_c_pbe(gga_type *p, double *rho, double *grho,
	       double *e, double *dedd, double *dedgd);


/* meta-GGAs */
void mgga_x_tpss_init(mgga_type *p);
void mgga_c_tpss_init(mgga_type *p);

void mgga_x_tpss_end(mgga_type *p);
void mgga_c_tpss_end(mgga_type *p);

void mgga_x_tpss(mgga_type *p, double *rho, double *grho, double *tau,
		 double *e, double *dedd, double *dedgd, double *dedtau);
void mgga_c_tpss(mgga_type *p, double *rho, double *grho, double *tau,
		 double *e, double *dedd, double *dedgd, double *dedtau);

/* LCAs */
void lca_lch_init(lca_type *p);
void lca_omc_init(lca_type *p);

void lca_s_lch(double rs, double s, double dsdrs);
void lca_s_omc(double rs, double s, double dsdrs);


#endif
