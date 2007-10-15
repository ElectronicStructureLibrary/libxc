#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_LCA_OMC       301   /* Orestes, Marcasso & Capelle  */

static xc_func_info_type func_info_lca_omc = {
  XC_LCA_OMC,
  XC_EXCHANGE_CORRELATION,
  "Orestes, Marcasso & Capelle parametrization",
  XC_FAMILY_LCA,
  "E. Orestes, T. Marcasso, and K. Capelle, Phys. Rev. A 68, 022105 (2003)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC
};

void lca_omc_init(xc_lca_type *p)
{
  p->info = &func_info_lca_omc;
}


/* This routine computes the ratio of the orbital susceptibilities for the interacting and 
   non-interacting electron gas and its derivative */
void lca_s_omc(double rs, double *s, double *dsdrs)
{
  static double c[5] = {1.1038, -0.4990, 0.4423, -0.06696, 0.0008432};
  double tmp;
  
  tmp    = sqrt(rs);
  *s     = c[0] + c[1]*pow(rs, 1.0/3.0) + c[2]*tmp + c[3]*rs + c[4]*rs*rs;
  *dsdrs = c[1]*pow(rs, -2.0/3.0)/3.0 + c[2]/(2.0*tmp) + c[3] + 2.0*c[4]*rs;
}
