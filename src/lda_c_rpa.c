#include <stdio.h>
#include "util.h"

/************************************************************************
 Random Phase Approximation (RPA)
************************************************************************/

#define XC_LDA_C_RPA  3   /* Random Phase Approximation   */

static void lda_c_rpa(const void *p_, const double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a = 0.0311, b = -0.047, c = 0.009, d = -0.017;
  double dens, zeta, rs;
  double lrs;

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs  =  RS(dens); /* Wigner radius */
  lrs = log(rs);
  
  *ec   = a*lrs + b + c*rs*lrs + d*rs;
  vc[0] = a/rs + c*(lrs + 1.0) + d;         /* now contains d ec/d rs */
  
  vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
  if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}

const xc_func_info_type func_info_lda_c_rpa = {
  XC_LDA_C_RPA,
  XC_CORRELATION,
  "Random Phase Approximation (RPA)",
  XC_FAMILY_LDA,
  "M. Gell-Mann and K.A. Brueckner, Phys. Rev. 106, 364 (1957)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_rpa,    /* lda  */
};
