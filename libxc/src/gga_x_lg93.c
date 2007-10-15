#include <stdio.h>
#include "util.h"

#define XC_GGA_X_LG93  113 /* Lacks & Gordon 93 */

static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  static const double ad = 1e-8, a4 = 29.790, a6 = 22.417;
  static const double a8 = 12.119, a10 = 1570.1, a12 = 55.944;
  static const double a2 = 4.94113918475214219939; /* (ad + 0.1234)/b, b = 0.024974 */

  double ss, ss2, ss4, ss6, ss8, ss10;
  double f1, f2, f3;

  ss  = x2s*x;    ss2  = ss*ss;
  ss4 = ss2*ss2;  ss6  = ss4*ss2;
  ss8 = ss6*ss2;  ss10 = ss8*ss2;

  f1 = 1.0 + a2*ss2 + a4*ss4 + a6*ss6 + a8*ss8 + a10*ss10 + a12*ss2*ss10;
  f2 = 1.0 + ad*ss2;

  *f = f1/f2;

  f3 = 2.0*ss*(a2 + 2.0*a4*ss2 + 3.0*a6*ss4 + 4.0*a8*ss6 + 5.0*a10*ss8 + 6.0*a12*ss10);
  *dfdx  = x2s*(f3*f2 - 2.0*ss*ad*f1)/(f2*f2);
  *ldfdx = x2s*x2s*(a2 - ad);
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_lg93 = {
  XC_GGA_X_LG93,
  XC_EXCHANGE,
  "Lacks & Gordon 93",
  XC_FAMILY_GGA,
  "D.J. Lacks and R.G. Gordon, Phys. Rev. A 47, 4681-4690 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

