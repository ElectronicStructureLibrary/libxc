#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Exchange part
************************************************************************/

static xc_func_info_type func_info_mgga_x_tpss = {
  XC_MGGA_X_TPSS,
  XC_EXCHANGE,
  "Perdew, Tao, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  "J.P.Perdew, Tao, Staroverov, and Scuseria, Phys. Rev. Lett. 91, 146401 (2003)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC
};


void mgga_x_tpss_init(xc_mgga_type *p)
{
  p->info = &func_info_mgga_x_tpss;
  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_x_init(p->lda_aux, XC_UNPOLARIZED, 3, XC_NON_RELATIVISTIC);
}


void mgga_x_tpss_end(xc_mgga_type *p)
{
  free(p->lda_aux);
}


/* some parameters */
static double b=0.40, c=1.59096, e=1.537, kappa=0.804, mu=0.21951;


/* This is Equation (7) from the paper and its derivatives */
static void 
x_tpss_7(double p, double z, 
	 double *qb, double *dqbdp, double *dqbdz)
{
  double alpha, dalphadp, dalphadz;

  { /* Eq. (8) */
    double a = (1.0/z - 1.0), h = 5.0/3.0;
    alpha    = h*a*p;
    dalphadp = h*a;
    dalphadz = -h*p/(z*z);
  }

  { /* Eq. (7) */
    double dqbda;
    double a = sqrt(1.0 + b*alpha*(alpha-1.0)), h = 9.0/20.0;
    dqbda = h*(1.0 + 0.5*b*(alpha-1.0))/pow(a, 3);

    *qb    = h*(alpha - 1.0)/a + 2.0*p/3.0;
    *dqbdp = dqbda*dalphadp + 2.0/3.0;
    *dqbdz = dqbda*dalphadz;
  }

}

/* Equation (10) in all it's glory */
static 
void x_tpss_10(double p, double z, 
	       double *x, double *dxdp, double *dxdz)
{
  double x1, dxdp1, dxdz1;
  double aux1, z2, p2;
  double qb, dqbdp, dqbdz;
  
  /* Equation 7 */
  x_tpss_7(p, z, &qb, &dqbdp, &dqbdz);

  z2   = z*z;
  p2   = p*p; 
  aux1 = 10.0/81.0;
  
  /* first we handle the numerator */
  x1    = 0.0;
  dxdp1 = 0.0;
  dxdz1 = 0.0;

  { /* first term */
    double a = 1.0+z2, a2 = a*a;
    x1    += (aux1 + c*z2/a2)*p;
    dxdp1 += (aux1 + c*z2/a2);
    dxdz1 += c*2.0*z*(1.0-z2)*p/(a*a2);
  }
  
  { /* second term */
    double a = 146.0/2025.0*qb;
    x1    += a*qb;
    dxdp1 += 2.0*a*dqbdp;
    dxdz1 += 2.0*a*dqbdz;
  }
  
  { /* third term */
    double a = sqrt(0.5*(9.0*z2/25.0 + p2));
    double h = 72.0/405;
    x1    += -h*qb*a;
    dxdp1 += -h*(a*dqbdp + 0.5*qb*p/a);
    dxdz1 += -h*(a*dqbdz + 0.5*qb*(9.0/25.0)*z/a);
  }
  
  { /* forth term */
    double a = aux1*aux1/kappa;
    x1    += a*p2;
    dxdp1 += a*2.0*p;
  }
  
  { /* fifth term */
    double a = 2.0*sqrt(e)*aux1*9.0/25.0;
    x1    += a*z2;
    dxdz1 += a*2.0*z;
  }
  
  { /* sixth term */
    double a = e*mu;
    x1    += a*p*p2;
    dxdp1 += a*3.0*p2;
  }
  
  /* and now the denominator */
  {
    double a = 1.0+sqrt(e)*p, a2 = a*a;
    *x    = x1/a2;
    *dxdp = (dxdp1*a - 2.0*sqrt(e)*x1)/(a2*a);
    *dxdz = dxdz1/a2;
  }
}

static void 
x_tpss_para(xc_mgga_type *pt, double rho, double *grho, double tau_,
	    double *energy, double *dedd, double *dedgd, double *dedtau)
{

  double gdms, p, tau, tauw, z;
  double x, dxdp, dxdz, Fx, dFxdx;
  double exunif, vxunif;
  
  tau = max(tau_, MIN_TAU);

  /* get the uniform gas energy and potential */
  xc_lda(pt->lda_aux, &rho, &exunif, &vxunif, NULL);

  /* calculate |nabla rho|^2 */
  gdms = grho[0]*grho[0] + grho[1]*grho[1] + grho[2]*grho[2];
  gdms = max(MIN_GRAD*MIN_GRAD, gdms);
  
  /* Eq. (4) */
  p = gdms/(4.0*pow(3*M_PI*M_PI, 2.0/3.0)*pow(rho, 8.0/3.0));

  /* von Weisaecker kinetic energy density */
  tauw = gdms/(8.0*rho);
  z  = tauw/tau;

  /* Eq. 10 */
  x_tpss_10(p, z, &x, &dxdp, &dxdz);

  { /* Eq. (5) */
    double a = kappa/(kappa + x);
    Fx    = 1.0 + kappa*(1.0 - a);
    dFxdx = a*a;
  }
  
  { /* Eq. (3) */
    int i;
    double a = rho*exunif*dFxdx;

    *energy = exunif*Fx;
    *dedd   = vxunif*Fx + exunif*dFxdx*(-(8.0/3.0)*p*dxdp - z*dxdz);
    *dedtau = a * (-z/tau*dxdz);

    for(i=0; i<3; i++)
      dedgd[i] = a * 2.0*grho[i]/gdms * (p*dxdp + z*dxdz);
  }
}


void 
mgga_x_tpss(xc_mgga_type *p, double *rho, double *grho, double *tau,
	    double *e, double *dedd, double *dedgd, double *dedtau)
{
  if(p->nspin == XC_UNPOLARIZED){
    x_tpss_para(p, rho[0], grho, tau[0], e, dedd, dedgd, dedtau);

  }else{ 
    /* The spin polarized version is handle using the exact spin scaling
          Ex[n1, n2] = (Ex[2*n1] + Ex[2*n2])/2
    */
    int is;

    *e = 0.0;
    for(is=0; is<2; is++){
      double gr[3], e1;
      int i;
      for(i=0; i<3; i++) gr[i] = 2.0*grho _(is, i);

      x_tpss_para(p, 2.0*rho[is], gr, 2.0*tau[is], &e1, 
		  &(dedd[is]), &(dedgd _(is, 0)), &(dedtau[is]));
      *e += e1;
    }
  }
}
