#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"


/************************************************************************
 Slater's Xalpha functional

    Exc = alpha Ex
************************************************************************/

/* This correlation functional, added to the exchange functional, produces
a total exchange and correlation functional, Exc, equal to 3/2 * alpha * Ex 
Setting alpha equal to one gives the *usual* Slater Xalpha functional,
whereas alpha equal to 2/3 just leaves the exchange functional unchanged */

#define XC_LDA_C_XALPHA  6   /* Slater's Xalpha              */

static void lda_c_xalpha(const void *p_, const double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;
  double a = 1.5*p->alpha - 1.0;
  int i;

  xc_lda(p->lda_aux, rho, ec, vc, fc, NULL);

  if(ec != NULL)
    (*ec) *= a;

  if(vc != NULL)
    for(i=0; i<p->nspin; i++) vc[i] *= a;

  if(fc != NULL){
    int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++) fc[i] *= a;
  }
}

/* These prototypes are needed for the declaration of func_info_lda_c_xalpha */
void xc_lda_c_xalpha_init_default(void *p_);
void xc_lda_c_xalpha_end(void *p_);

const xc_func_info_type func_info_lda_c_xalpha = {
  XC_LDA_C_XALPHA,
  XC_CORRELATION,
  "Slater's Xalpha",
  XC_FAMILY_LDA,
  NULL,
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  xc_lda_c_xalpha_init_default,  /* init */
  xc_lda_c_xalpha_end,           /* end  */
  lda_c_xalpha                   /* lda */
};

void xc_lda_c_xalpha_init(xc_lda_type *p, int nspin, int dim, double alpha)
{
  p->info = &func_info_lda_c_xalpha;
  p->nspin = nspin;
  p->dim   = dim;
  p->alpha = alpha;

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_x_init(p->lda_aux, nspin, dim, XC_NON_RELATIVISTIC);
}

void xc_lda_c_xalpha_init_default(void *p_)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  xc_lda_c_xalpha_init(p, p->nspin, 3, 1.0);
}

void xc_lda_c_xalpha_end(void *p_)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  free(p->lda_aux);
}
