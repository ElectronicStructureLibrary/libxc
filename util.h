#ifndef _LDA_H
#define _LDA_H

#include <math.h>
#include "xc.h"

#define min(x,y)  ((x<y) ? x : y)
#define max(x,y)  ((x<y) ? y : x)
#define RS(x)     (pow((3.0/(4*M_PI*x)), 1.0/3.0))
#define FZETA(x)  ((pow(1.0 + zeta, 4.0/3.0) + pow(1.0 - zeta, 4.0/3.0) - 2.0)/0.519842099789746380)
#define DFZETA(x) ((pow(1.0 + zeta, 1.0/3.0) - pow(1.0 - zeta, 1.0/3.0))*(4.0/3.0)/0.519842099789746380)

#define MIN_DENS             1.0e-14
#define MIN_GRAD             1.0e-14

#define _(is, x) [3*is + x]

void rho2dzeta(int nspin, double *rho, double *d, double *zeta);

/* LDAs */
void lda_x       (lda_type *p, double *rho, double *ex, double *vx);
void lda_c_wigner(lda_type *p, double rs, double *ec, double *vc);
void lda_c_rpa   (lda_type *p, double rs, double *ec, double *vc);
void lda_c_hl    (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_xalpha(lda_type *p, double *rho, double *ec, double *vc);
void lda_c_vwn_init();
void lda_c_vwn   (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_pz    (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_pw    (lda_type *p, double rs, double zeta, double *ec, double *vc);
void lda_c_amgb_init();
void lda_c_amgb  (lda_type *p, double *rho, double *ec, double *vc);


/* GGAs */
void gga_x_pbe_init(gga_type *p);
void gga_c_pbe_init(gga_type *p);
void gga_x_pbe     (gga_type *p, double *rho, double *grho,
		    double *e, double *dedd, double *dedgd);
void gga_c_pbe     (gga_type *p, double *rho, double *grho,
		    double *e, double *dedd, double *dedgd);
void gga_pbe_end (gga_type *p);
#endif
