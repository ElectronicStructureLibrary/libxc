#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"


static func_type func_lca_lch = {
  XC_LCA_LCH,
  XC_EXCHANGE_CORRELATION,
  "Lee, Colwell & Handy parametrization",
  "Vignale-Rasolt CDFT functional",
  "A.M. Lee, S.M. Colwell and N.C. Handy, Chem. Phys. Lett. 229, 225 (1994)"
  "A.M. Lee, S.M. Colwell and N.C. Handy, J. Chem. Phys. 103, 10095 (1995)"
};


void lca_lch_init(lca_type *p)
{
  p->func = &func_lca_lch;
}


/* This routine computes the ratio of the orbital susceptibilities fo the interactiog and 
   non-interacting electron gas and its derivative */
void lca_s_lch(double rs, double s, double dsdrs)
{
  static double c[3] = {1.0, 0.028, -0.042};
  double tmp;

  tmp   = exp(c[2]*rs);
  s     = (c[0] + c[1]*rs)*tmp;
  dsdrs = (c[1] + c[2]*(c[0] + c[1]*rs))*tmp;
}
