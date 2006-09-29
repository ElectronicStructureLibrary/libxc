#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Becke 86 Generalized Gradient Approximation.
************************************************************************/

static void gga_x_b86(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  gga_type *p = p_;

  static const double c     = 0.9305257363491; /* 3/8*cur(3/pi)*4^(2/3) */
  static const double beta[2]  = {
    0.0076,  /* beta from the original Becke paper */
    0.00787  /* reoptimized value used in part 3 os Becke 1997 paper */
  };
  static const double gamma = 0.004;

  double sfact, dens;
  int is;
  int func = p->func->number - XC_GGA_X_B86;

  assert(func==0 || func==1);

  *e   = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact     = 1.0;
    vsigma[1] = 0.0; /* there are no cross terms in this functional */
  }else
    sfact     = 2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    double gdm;
    double x, f1, f2, f, ds, rho13;
    double dfdx;
    int js = is==0 ? 0 : 2;

    if(rho[is] < MIN_DENS){
      vrho[is] = 0.0;
      vsigma[js] = 0.0;
      continue;
    }

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = pow(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);
    f1    = (1.0 + beta[func]*x*x);
    f2    = (1.0 + gamma*x*x);
    f     = f1/f2;
    
    (*e) += -sfact*c*(ds*rho13)*f;
      
    dfdx = 2.0*x*(beta[func]*f2 - gamma*f1)/(f2*f2);
      
    vrho[is]   = -4.0/3.0*c*rho13*(f - dfdx*x);
    if(gdm>MIN_GRAD)
      vsigma[js] = -sfact*c*(ds*rho13)*dfdx*x/(2.0*sigma[js]);
    else /* this is the correct limit, I think */
      vsigma[js] = -c/(sfact*(ds*rho13))*(beta[func]-gamma);
  }

  *e /= dens; /* we want energy per particle */
}

func_type func_gga_x_b86 = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  "GGA",
  "A. D. Becke, J. Chem. Phys 84, 4524 (1986)",
  NULL,
  NULL,
  NULL,
  gga_x_b86
};

func_type func_gga_x_b86_r = {
  XC_GGA_X_B86_R,
  XC_EXCHANGE,
  "Becke 86 (reoptimized)",
  "GGA",
  "A. D. Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "A. D. Becke, J. Chem. Phys 107, 8554 (1997)",
  NULL,
  NULL,
  NULL,
  gga_x_b86
};
