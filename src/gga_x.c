#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************/

static void pbe_f(int func, double x, double *f, double *dfdx, double *ldfdx)
{
  static const  double kappa[2] = {
    0.8040,
    1.245
  };
  static const double mu = 0.00361218645365094697; /* beta*(pi^2/(3^5*2^8))^(1/3) */
  double dd;

  dd     = 1.0/(kappa[func] + mu*x*x);

  *f     = 1.0 + kappa[func]*(1.0 - kappa[func]*dd);
  *dfdx  = 2.0*x*mu*kappa[func]*kappa[func]*dd*dd;
  *ldfdx = mu;
}


static void pw86_f(double x, double *f, double *dfdx, double *ldfdx)
{
  const double x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  const double aa = 1.296, bb = 14.0, cc = 0.2;
  double ss, ss2, ss4, dd;

  ss  = x2s*x;
  ss2 = ss*ss;
  ss4 = ss2*ss2;

  dd = 1.0 + aa*ss2 + bb*ss4 + cc*ss4*ss2;

  *f     = pow(dd, 1.0/15.0);
  *dfdx  = x2s*ss*(2.0*aa + 4.0*bb*ss2 + 6.0*cc*ss4)/15.0 * pow(dd, -14.0/15.0);
  *ldfdx = x2s*x2s*aa/15.0;
}


static void pw91_f(double x, double *f, double *dfdx, double *ldfdx)
{
  const double x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  const double aa = 0.19645, bb = 7.7956, cc = 0.2743, dd=-0.1508, ff=0.004, alpha=100.0;
  double ss, ss2, ss4;
  double f1, f2, f3, f4;

  ss  = x2s*x;
  ss2 = ss*ss;
  ss4 = ss2*ss2;

  f1 = dd*exp(-alpha*ss2);
  f2 = aa*asinh(bb*ss);
  f3 = (cc + f1)*ss2 - ff*ss4;
  f4 = 1.0 + ss*f2 + ff*ss4;

  *f     = 1.0 + f3/f4;
  *dfdx  = (2.0*ss*(cc + f1*(1.0 - alpha*ss2) - 2.0*ff*ss2)*f4 - 
	    f3*(f2 + ss*(aa*bb/sqrt(1.0 + bb*bb*ss2) + 4.0*ff*ss2)))/(f4*f4);
  *dfdx *= x2s;
  *ldfdx = x2s*x2s*(cc + dd);
}


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


static void optx_f(double x, double *f, double *dfdx, double *ldfdx)
{
  static const double a = 1.05151, b = 1.43169/X_FACTOR_C, gamma = 0.006;

  double f1, u;

  f1 = 1.0 + gamma*x*x;
  u  = gamma*x*x/f1;

  *f     = a + b*u*u;
  *dfdx  = 2.0*b*u * 2.0*gamma*x/(f1*f1);
  *ldfdx = 0.0;
}


static void dk87_f(int func, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double a1[2] = {0.861504, 0.861213}, 
    b1[2] = {0.044286, 0.042076}, alpha[2] = {1.0, 0.98};
  static const double betag = 0.00132326681668994855/X_FACTOR_C; /* 7/(432*pi*(6*pi^2)^(1/3)) */
  
  double f0, f1, f2;

  f0 = a1[func]*pow(x, alpha[func]);
  f1 = 1.0 + f0;
  f2 = 1.0 + b1[func]*x*x;
  
  *f     = 1.0 + betag*x*x*f1/f2;
  *dfdx  = betag*(2.0*x*f1/f2 + x*(alpha[func]*f0*f2 - 2.0*b1[func]*x*x*f1)/(f2*f2));
  *ldfdx = betag;
}


static void lg93_f(double x, double *f, double *dfdx, double *ldfdx)
{
  static const double x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  static const double ad = 1e-8, b = 0.024974, a4 = 29.790, a6 = 22.417;
  static const double a8 = 12.119, a10 = 1570.1, a12 = 55.944;
  static const double a2 = 4.94113918475214219939; /* (ad + 0.1234)/b */

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


/************************************************************************/

void gga_x_b86(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  xc_gga_type *p = p_;

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

    switch(p->info->number){
    case XC_GGA_X_PBE:
      pbe_f(0, x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_PBE_R:
      pbe_f(1, x, &f, &dfdx, &ldfdx);
      break;
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
    case XC_GGA_X_PW86:
      pw86_f(x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_PW91:
      pw91_f(x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_OPTX:
      optx_f(x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_DK87_R1:
      dk87_f(0, x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_DK87_R2:
      dk87_f(1, x, &f, &dfdx, &ldfdx);
      break;
    case XC_GGA_X_LG93:
      lg93_f(x, &f, &dfdx, &ldfdx);
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

const xc_func_info_type func_info_gga_x_pbe = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_pbe_r = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)\n"
  "Y. Zhang and W. Yang, Phys. Rev. Lett 80, 890 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_b86 = {
  XC_GGA_X_B86,
  XC_EXCHANGE,
  "Becke 86",
  XC_FAMILY_GGA,
  "A.D. Becke, J. Chem. Phys 84, 4524 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
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
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_b86_mgc = {
  XC_GGA_X_B86_MGC,
  XC_EXCHANGE,
  "Becke 86 with modified gradient correction",
  XC_FAMILY_GGA,
  "A.D. Becke, J. Chem. Phys 84, 4524 (1986)\n"
  "A.D. Becke, J. Chem. Phys 85, 7184 (1986)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  "A.D. Becke, Phys. Rev. A 38, 3098-3100 (1988)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_g96 = {
  XC_GGA_X_G96,
  XC_EXCHANGE,
  "Gill 96",
  XC_FAMILY_GGA,
  "P.M.W. Gill, Mol. Phys. 89, 433 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_pw86 = {
  XC_GGA_X_PW86,
  XC_EXCHANGE,
  "Perdew & Wang 86",
  XC_FAMILY_GGA,
  "J.P. Perdew and Y. Wang, Phys. Rev. B 33, 8800 (1986).",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_pw91 = {
  XC_GGA_X_PW91,
  XC_EXCHANGE,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  "J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, M.R. Pederson, and C. Fiolhais, Phys. Rev. B 46, 6671 (1992)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_optx = {
  XC_GGA_X_OPTX,
  XC_EXCHANGE,
  "Handy & Cohen OPTX 01",
  XC_FAMILY_GGA,
  "N.C. handy and A.J. Cohen, Mol. Phys. 99, 403-412 (2001)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_dk87_r1 = {
  XC_GGA_X_DK87_R1,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R1",
  XC_FAMILY_GGA,
  "A.E. DePristo and J.D. Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_dk87_r2 = {
  XC_GGA_X_DK87_R2,
  XC_EXCHANGE,
  "dePristo & Kress 87 version R2",
  XC_FAMILY_GGA,
  "A.E. DePristo and J.D. Kress, J. Chem. Phys. 86, 1425 (1987)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};

const xc_func_info_type func_info_gga_x_lg93 = {
  XC_GGA_X_LG93,
  XC_EXCHANGE,
  "Lacks & Gordon 93",
  XC_FAMILY_GGA,
  "D.J. Lacks and R.G. Gordon, Phys. Rev. A 47, 4681-4690 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_x_b86
};
