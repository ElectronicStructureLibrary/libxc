#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Perdew, Burke & Ernzerhof Generalized Gradient Approximation.

   [1] J.P.Perdew, K.Burke & M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)

 I based this implementation on a rotine from L.C. Balbas and J.M. Soler
************************************************************************/

/* some parameters */
static const double beta = 0.066725, kappa = 0.8040;
static const double gamm = 0.0310906908696549; /* (1.0 - log(2.0))/(M_PI*M_PI) */
static const double mu   = 0.2195164512208958; /* beta*M_PI*M_PI/3.0 */

void gga_x_pbe_init(gga_type *p)
{
  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_x_init(p->lda_aux, XC_UNPOLARIZED, 3, XC_NON_RELATIVISTIC);
}


void gga_c_pbe_init(gga_type *p)
{
  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_init(p->lda_aux, XC_LDA_C_PW, p->nspin);
}


void gga_pbe_end(gga_type *p)
{
  free(p->lda_aux);
}


void gga_x_pbe(gga_type *p, double *rho, double *grho,
	       double *e, double *dedd, double *dedgd)
{
  double dens, sfact;
  int is, ix;

  *e   = 0.0;
  dens = 0.0;
  sfact = (p->nspin == XC_POLARIZED) ? 2.0 : 1.0;

  for(is=0; is<p->nspin; is++){
    double ds, gdms, kfs, s, f1, f;
    double exunif, vxunif;
    double dkfdd, dsdd, df1dd, dfdd;

    dens = dens + rho[is];
    ds   = max(MIN_DENS, sfact*rho[is]);

    /* calculate |nabla rho| */
    gdms = sqrt(grho _(is, 0)*grho _(is, 0) + 
		grho _(is, 1)*grho _(is, 1) +
		grho _(is, 2)*grho _(is, 2));
    gdms = max(MIN_GRAD, sfact*gdms);

    kfs  = pow(3.0*M_PI*M_PI*ds, 1.0/3.0);
    s    = gdms/(2.0*kfs*ds);

    f1   = 1.0 + mu*s*s/kappa;
    f    = 1.0 + kappa - kappa/f1;

    lda(p->lda_aux, &ds, &exunif, &vxunif);
    
    /* total energy per unit volume */
    *e += ds*exunif*f;
    
    dkfdd = kfs/(3.0*ds);
    dsdd  = s*(-dkfdd/kfs - 1.0/ds);
    df1dd = 2.0*(f1 - 1.0)*dsdd/s;
    dfdd  = kappa*df1dd/(f1*f1);

    dedd[is] = vxunif*f + ds*exunif*dfdd;
    
    for(ix=0; ix<3; ix++){
      double gds, dsdgd, df1dgd, dfdgd;

      gds    = sfact*grho _(is, ix);
      dsdgd  = (s/gdms)*gds/gdms;
      df1dgd = 2.0*mu*s*dsdgd/kappa;
      dfdgd  = kappa*df1dgd/(f1*f1);

      dedgd _(is,ix) = ds*exunif*dfdgd;
    }
  }

  *e = *e/(dens*sfact); /* we want energy per particle */
}

void gga_c_pbe(gga_type *p, double *rho, double *grho,
	       double *e, double *dedd, double *dedgd)
{
  double dens, zeta, ecunif, vcunif[2];
  double rs, kf, ks, phi, phi3, gdmt, t, t2;
  double f1, f2, f3, f4, a, h;
  double drsdd, dkfdd, dksdd, dzdd[2], dpdz;
  int is, ix;

  lda(p->lda_aux, rho, &ecunif, vcunif);
  rho2dzeta(p->nspin, rho, &dens, &zeta);
  
  rs = RS(dens);
  kf = pow(3.0*M_PI*M_PI*dens, 1.0/3.0);
  ks = sqrt(4.0*kf/M_PI);

  phi  = 0.5*(pow(1.0 + zeta, 2.0/3.0) + pow(1.0 - zeta, 2.0/3.0));
  phi3 = pow(phi, 3);

  /* get gdmt = |nabla n| */
  gdmt = 0.0;
  for(ix=0; ix<3; ix++){
    double dd = grho _(0, ix);
    if(p->nspin == XC_POLARIZED) dd += grho _(1, ix);
    gdmt += dd*dd;
  }
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
  dzdd[1] =  (1.0 - zeta)/dens;
  dzdd[2] = -(1.0 + zeta)/dens;
  dpdz    = (1.0/3.0)*(1/pow(1 + zeta, 1.0/3.0) - 1/pow(1 - zeta, 1.0/3.0));
  
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
    dedd[is] = vcunif[is] + h + dens*dhdd;

    for(ix=0; ix<3; ix++){
      double gdt, dtdgd, df3dgd, df4dgd, dhdgd;

      gdt = grho _(0, ix);
      if(p->nspin == XC_POLARIZED) gdt += grho _(1, ix);

      dtdgd  = t*gdt/(gdmt*gdmt);
      df3dgd = dtdgd*t*(2.0 + 4.0*a*t2);
      df4dgd = f4*df3dgd*(1.0/f3 - a/(1.0 + a*f3));
      dhdgd  = gamm*phi3*df4dgd/(1.0 + f4);
      dedgd _(is, ix) = dens*dhdgd;
    }
  }
   
}
