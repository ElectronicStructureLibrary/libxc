#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************/

static void b86_f(int func, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double beta[2]  = {
    0.0076,  /* beta from the original Becke paper */
    0.00787  /* reoptimized value used in part 3 of Becke 97 paper */
  };
  static const double gamma = 0.004;

  double f1, f2;

  f1    = (1.0 + beta[func]*x*x);
  f2    = (1.0 + gamma*x*x);
  *f    = f1/f2;
    
  *dfdx  = 2.0*x*(beta[func]*f2 - gamma*f1)/(f2*f2);
  *ldfdx = (beta[func] - gamma);
}


static void b86_mgc_f(double x, double *f, double *dfdx, double *ldfdx)
{
  static const double beta  = 0.00375;
  static const double gamma = 0.007;
  
  double f1;

  f1    = (1.0 + gamma*x*x);
  *f    = 1.0 + beta/X_FACTOR_C*x*x/pow(f1, 4.0/5.0);

  *dfdx = beta/X_FACTOR_C*2.0*x*(5.0 + gamma*x*x)/(5.0*pow(f1, 9.0/5.0));
  *ldfdx= beta/X_FACTOR_C;

}


static void b88_f(double x, double *f, double *dfdx, double *ldfdx)
{
  static const double beta  = 0.0042;

  double f1;

  f1 = (1.0 + 6.0*beta*x*asinh(x));
  *f = 1.0 + beta/X_FACTOR_C*x*x/f1;
 
  *dfdx = beta/X_FACTOR_C*x*(2.0 + 6.0*beta*(x*asinh(x) - x*x/sqrt(1.0+x*x)))/(f1*f1);
  *ldfdx= beta/X_FACTOR_C;
}


static void g96_f(double x, double *f, double *dfdx, double *ldfdx)
{
  static const double c1 = 1.0/137.0;
  double sx = sqrt(x);

  *f     = 1.0 + c1/X_FACTOR_C*x*sx;
  *dfdx  = 3.0*c1/(2.0*X_FACTOR_C)*sx;
  *ldfdx = 0.0; /* This is not true, but I think this functional diverges */
}


/************************************************************************/

void gga_x_b86(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  gga_type *p = p_;

  double sfact, dens;
  int is;

  *e   = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact     = 1.0;
    vsigma[1] = 0.0; /* there are no cross terms in this functional */
  }else
    sfact     = 2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    double gdm, ds, rho13;
    double x, f, dfdx, ldfdx;
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

    switch(p->func->number){
    case XC_GGA_X_B86:
      b86_f(0, x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_B86_R:
      b86_f(1, x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_B86_MGC:
      b86_mgc_f(x, &f, &dfdx, &ldfdx);
      break; 
    case XC_GGA_X_B88:
      b88_f(x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_G96:
      g96_f(x, &f, &dfdx, &ldfdx);
      break;
   default:
      abort();
    }
    
    (*e) += -sfact*X_FACTOR_C*(ds*rho13)*f;
      
    vrho[is]   = -4.0/3.0*X_FACTOR_C*rho13*(f - dfdx*x);
    if(gdm>MIN_GRAD)
      vsigma[js] = -sfact*X_FACTOR_C*(ds*rho13)*dfdx*x/(2.0*sigma[js]);
    else
      vsigma[js] = -X_FACTOR_C/(sfact*(ds*rho13))*ldfdx;
  }

  *e /= dens; /* we want energy per particle */
}


/************************************************************************/

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


func_type func_gga_x_b86_mgc = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  "GGA",
  "A. D. Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "A. D. Becke, J. Chem. Phys 85, 7184 (1986)",
  NULL,
  NULL,
  NULL,
  gga_x_b86
};

const func_type func_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  "GGA",
  "A. D. Becke, Phys. Rev. A 38, 3098-3100 (1988)",
  NULL,
  NULL,
  NULL,
  gga_x_b86
};

const func_type func_gga_x_g96 = {
  XC_GGA_X_G96,
  XC_EXCHANGE,
  "Gill 96",
  "GGA",
  "P. M. W. Gill, Mol. Phys. 89, 433 (1996)",
  NULL,
  NULL,
  NULL,
  gga_x_b86
};
