#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_RPBE  117 /* Hammer, Hansen & Norskov (PBE-like) */

/* RPBE: see PBE for more details */
static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double kappa = 0.8040;
  static const double mu = 0.00361218645365094697;

  double dd;

  dd     = exp(-mu*x*x/kappa);

  *f     = 1.0 + kappa*(1.0 - dd);
  *dfdx  = 2.0*x*mu*dd;
  *ldfdx = mu;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_rpbe = {
  XC_GGA_X_RPBE,
  XC_EXCHANGE,
  "Hammer, Hansen, and Nørskov",
  XC_FAMILY_GGA,
  "B Hammer, LB Hansen and JK Nørskov, Phys. Rev. B 59, 7413 (1999)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
