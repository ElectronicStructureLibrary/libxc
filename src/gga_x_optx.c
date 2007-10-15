#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_OPTX         110 /* Handy & Cohen OPTX 01                          */

static inline void
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double a = 1.05151, b = 1.43169/X_FACTOR_C, gamma = 0.006;

  double f1, u;

  f1 = 1.0 + gamma*x*x;
  u  = gamma*x*x/f1;

  *f     = a + b*u*u;
  *dfdx  = 2.0*b*u * 2.0*gamma*x/(f1*f1);
  *ldfdx = 0.0;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  "N.C. handy and A.J. Cohen, Mol. Phys. 99, 403-412 (2001)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
