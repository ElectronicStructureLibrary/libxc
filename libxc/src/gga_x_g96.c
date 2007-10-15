#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_G96          107 /* Gill 96                                        */

static inline void
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double c1 = 1.0/137.0;
  double sx = sqrt(x);

  *f     = 1.0 + c1/X_FACTOR_C*x*sx;
  *dfdx  = 3.0*c1/(2.0*X_FACTOR_C)*sx;
  *ldfdx = 0.0; /* This is not true, but I think this functional diverges */
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_g96 = {
  XC_GGA_X_G96,
  XC_EXCHANGE,
  "Gill 96",
  XC_FAMILY_GGA,
  "P.M.W. Gill, Mol. Phys. 89, 433 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
