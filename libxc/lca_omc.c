#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"


static func_type func_lca_omc = {
  XC_LCA_OMC,
  XC_EXCHANGE_CORRELATION,
  "Orestes, Marcasso & Capelle parametrization",
  "Vignale-Rasolt CDFT functional",
  "E. Orestes, T. Marcasso, and K. Capelle, Phys. Rev. A 68, 022105 (2003)"
};

void lca_omc_init(lca_type *p)
{
  p->func = &func_lca_omc;
}


/* This routine computes the ratio of the orbital susceptibilities fo the interactiog and 
   non-interacting electron gas and its derivative */
void lca_s_omc(double rs, double s, double dsdrs)
{
  static double c[5] = {1.1038, -0.4990, 0.4423, -0.06696, 0.0008432};
  double tmp;
  
  tmp   = sqrt(rs);
  s     = c[0] + c[1]*pow(rs, 1.0/3.0) + c[2]*tmp + c[3]*rs + c[4]*rs*rs;
  dsdrs = c[1]*pow(rs, -2.0/3.0)/3.0 + c[2]/(2.0*tmp) + c[3] + 2.0*c[4]*rs;
}
