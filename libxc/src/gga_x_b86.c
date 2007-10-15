#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_B86          103 /* Becke 86 Xalfa,beta,gamma                      */
#define XC_GGA_X_B86_R        104 /* Becke 86 Xalfa,beta,gamma (reoptimized)        */

static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double beta[2]  = {
    0.0076,  /* beta from the original Becke paper */
    0.00787  /* reoptimized value used in part 3 of Becke 97 paper */
  };
  static const double gamma = 0.004;

  double f1, f2;
  int func;

  switch(p->info->number){
  case XC_GGA_X_B86_R:   func = 1; break;
  default:               func = 0; /* original B86 */
  }

  f1    = (1.0 + beta[func]*x*x);
  f2    = (1.0 + gamma*x*x);
  *f    = f1/f2;
    
  *dfdx  = 2.0*x*(beta[func]*f2 - gamma*f1)/(f2*f2);
  *ldfdx = (beta[func] - gamma);
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_b86 = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  XC_FAMILY_GGA,
  "A.D. Becke, J. Chem. Phys 84, 4524 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

const xc_func_info_type func_info_gga_x_b86_r = {
  XC_GGA_X_B86_R,
  XC_EXCHANGE,
  "Becke 86 (reoptimized)",
  XC_FAMILY_GGA,
  "A.D. Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "A.D. Becke, J. Chem. Phys 107, 8554 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

