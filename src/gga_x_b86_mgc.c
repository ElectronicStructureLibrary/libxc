#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalfa,beta,gamma (with mod. grad. correction) */

static inline void
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double beta  = 0.00375;
  static const double gamma = 0.007;
  
  double f1;

  f1    = (1.0 + gamma*x*x);
  *f    = 1.0 + beta/X_FACTOR_C*x*x/pow(f1, 4.0/5.0);

  *dfdx = beta/X_FACTOR_C*2.0*x*(5.0 + gamma*x*x)/(5.0*pow(f1, 9.0/5.0));
  *ldfdx= beta/X_FACTOR_C;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_b86_mgc = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  "A.D. Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "A.D. Becke, J. Chem. Phys 85, 7184 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
