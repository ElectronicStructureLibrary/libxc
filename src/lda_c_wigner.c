#include <stdio.h>
#include "util.h"

/************************************************************************
 Wigner's parametrization from the low density limit
************************************************************************/

#define XC_LDA_C_WIGNER  2   /* Wigner parametrization       */

static void lda_c_wigner(const void *p_, const double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a = -0.44, b = 7.8;
  double dens, zeta, rs;
  double etmp, decdrs, t;
  
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs    =  RS(dens); /* Wigner radius */
  t     =  b + rs;

  etmp   =  a/t;
  decdrs = -a/(t*t);                         /* now contains d ec/d rs */
  
  if(ec != NULL) *ec = etmp;

  if(vc != NULL){
    vc[0] = etmp - decdrs*rs/3.0;              /* and now d ec/d rho */
    if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to return something */
  }

}

const xc_func_info_type func_info_lda_c_wigner = {
  XC_LDA_C_WIGNER,
  XC_CORRELATION,
  "Wigner",
  XC_FAMILY_LDA,
  "E.P. Wigner, Trans. Faraday Soc. 34, 678 (1938)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_wigner, /* lda  */
};
