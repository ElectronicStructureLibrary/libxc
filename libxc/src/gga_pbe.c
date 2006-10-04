#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Perdew, Burke & Ernzerhof Generalized Gradient Approximation.

 I based this implementation on a routine from L.C. Balbas and J.M. Soler
************************************************************************/


/*                            exchange                                 */
static const double kappa[2] = {
  0.8040,
  1.245
};
static const double mu    = 0.2195149727645171;  /* beta*M_PI*M_PI/3.0 */

void gga_x_pbe_init(void *p_)
{
  gga_type *p = p_;

  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_x_init(p->lda_aux, XC_UNPOLARIZED, 3, XC_NON_RELATIVISTIC);
}


void gga_x_pbe_end(void *p_)
{
  gga_type *p = p_;

  free(p->lda_aux);
}


void gga_x_pbe(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  gga_type *p = p_;

  double dens, sfact;
  int is;
  int func = p->func->number - XC_GGA_X_PBE;
  
  assert(func==0 || func==1);

  *e   = 0.0;
  dens = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact     = 2.0;
    vsigma[1] = 0.0; /* there are no cross terms in this functional */
  }else
    sfact     = 1.0;

  for(is=0; is<p->nspin; is++){
    double ds, gdms, kfs, s, f1, f, sig;
    double exunif, vxunif;
    double dkfdd, dsdd, df1dd, dfdd;
    double dsdsig, df1dsig, dfdsig;

    dens = dens + rho[is];
    ds   = max(MIN_DENS, sfact*rho[is]);

    /* calculate |nabla rho| */
    sig  = sfact*sfact*sigma[is==0 ? 0 : 2];
    gdms = max(MIN_GRAD, sqrt(sig));

    kfs  = pow(3.0*M_PI*M_PI*ds, 1.0/3.0);
    s    = gdms/(2.0*kfs*ds);

    f1   = 1.0 + mu*s*s/kappa[func];
    f    = 1.0 + kappa[func] - kappa[func]/f1;

    lda(p->lda_aux, &ds, &exunif, &vxunif, NULL);
    
    /* total energy per unit volume */
    *e += ds*exunif*f;
    
    dkfdd = kfs/(3.0*ds);
    dsdd  = s*(-dkfdd/kfs - 1.0/ds);
    df1dd = 2.0*(f1 - 1.0)*dsdd/s;
    dfdd  = kappa[func]*df1dd/(f1*f1);

    vrho[is] = vxunif*f + ds*exunif*dfdd;
    
    dsdsig  = s/(2.0*sig);
    df1dsig = 2.0*mu*s*dsdsig/kappa[func];
    dfdsig  = kappa[func]*df1dsig/(f1*f1);
    vsigma[is==0 ? 0 : 2] = sfact*ds*exunif*dfdsig;
  }

  *e = *e/(dens*sfact); /* we want energy per particle */
}


const func_type func_gga_x_pbe = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_x_pbe_init,
  gga_x_pbe_end,
  NULL,            /* this is not an LDA                   */
  gga_x_pbe,
};

const func_type func_gga_x_pbe_r = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "Y. Zhang and W. Yang, Phys. Rev. Lett 80, 890 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_x_pbe_init,
  gga_x_pbe_end,
  NULL,            /* this is not an LDA                   */
  gga_x_pbe,
};

/*                            correlation                                 */
static const double beta  = 0.06672455060314922;
static const double gamm  = 0.03109076908696549; /* (1.0 - log(2.0))/(M_PI*M_PI) */

void gga_c_pbe_init(void *p_)
{
  gga_type *p = p_;

  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_init(p->lda_aux, XC_LDA_C_PW, p->nspin);
}


void gga_c_pbe(void *p_, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  gga_type *p = p_;

  double dens, zeta, ecunif, vcunif[2];
  double rs, kf, ks, phi, phi3, gdmt, t, t2;
  double f1, f2, f3, f4, a, h;
  double drsdd, dkfdd, dksdd, dzdd[2], dpdz;
  int is;

  lda(p->lda_aux, rho, &ecunif, vcunif, NULL);
  rho2dzeta(p->nspin, rho, &dens, &zeta);
  
  rs = RS(dens);
  kf = pow(3.0*M_PI*M_PI*dens, 1.0/3.0);
  ks = sqrt(4.0*kf/M_PI);

  phi  = 0.5*(pow(1.0 + zeta, 2.0/3.0) + pow(1.0 - zeta, 2.0/3.0));
  phi3 = pow(phi, 3);

  /* get gdmt = |nabla n| */
  gdmt = sigma[0];
  if(p->nspin == XC_POLARIZED) gdmt += 2.0*sigma[1] + sigma[2];
  gdmt = sqrt(gdmt);

  t  = gdmt/(2.0*phi*ks*dens);
  t2 = t*t;

  f1 = ecunif/(gamm*phi3);
  f2 = exp(-f1);
  a  = beta/(gamm*(f2 - 1.0));
  f3 = t2 + a*t2*t2;
  f4 = beta*f3/(gamm*(1.0 + a*f3));
  h  = gamm*phi3*log(1.0 + f4);
  *e = ecunif + h;

  drsdd   = -rs/(3.0*dens);
  dkfdd   =  kf/(3.0*dens);
  dksdd   = 0.5*ks*dkfdd/kf;
  dzdd[0] =  (1.0 - zeta)/dens;
  dzdd[1] = -(1.0 + zeta)/dens;
  dpdz    = 0.0;
  if(fabs(1.0 + zeta) >= MIN_DENS)
    dpdz += (1.0/3.0)/pow(1.0 + zeta, 1.0/3.0);
  if(fabs(1.0 - zeta) >= MIN_DENS)
    dpdz -= (1.0/3.0)/pow(1.0 - zeta, 1.0/3.0);
  
  for(is=0; is<p->nspin; is++){
    double decudd, dpdd, dtdd;
    double df1dd, df2dd, df3dd, df4dd, dadd, dhdd;

    decudd = (vcunif[is] - ecunif)/dens;
    dpdd   = dpdz*dzdd[is];
    dtdd   = (-t)*(dpdd/phi + dksdd/ks + 1.0/dens);
    df1dd  = f1*(decudd/ecunif - 3.0*dpdd/phi);
    df2dd  = (-f2)*df1dd;
    dadd   = (-a)*df2dd/(f2 - 1.0);
    df3dd  = t*(2.0 + 4.0*a*t2)*dtdd + dadd*t2*t2;
    df4dd  = f4*(df3dd/f3 - (dadd*f3 + a*df3dd)/(1.0 + a*f3));
    dhdd   = 3.0*h*dpdd/phi;
    dhdd  += gamm*phi3*df4dd/(1.0 + f4);
    vrho[is] = vcunif[is] + h + dens*dhdd;
  }

  { /* calculate now vsigma */
    double dtdsig, df3dsig, df4dsig, dhdsig;
    
    dtdsig  = t/(2.0*gdmt*gdmt);
    df3dsig = dtdsig*t*(2.0 + 4.0*a*t2);
    df4dsig = f4*df3dsig*(1.0/f3 - a/(1.0 + a*f3));
    dhdsig  = gamm*phi3*df4dsig/(1.0 + f4);
    vsigma[0] = dens*dhdsig;
    if(is == 2){
      vsigma[1] = 2.0*vsigma[0];
      vsigma[2] =     vsigma[0];
    }
  }
}

const func_type func_gga_c_pbe = {
  XC_GGA_C_PBE,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_LDA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_c_pbe_init,
  gga_x_pbe_end,   /* we can use the same as exchange here */
  NULL,            /* this is not an LDA                   */
  gga_c_pbe,
};
