#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

void gga_lb_end(void *p_)
{
  gga_type *p = p_;

  free(p->lda_aux);
}

void gga_lb(gga_type *p, double *rho, double *sigma, double r, double ip, double qtot,
	    double *dedd)
{
  int is;
  double alpha, gdm, gamm, x;

  static const double beta = 0.05;

  lda(p->lda_aux, rho, &x, dedd, NULL);
  if(p->modified){
    alpha = (ip > 0.0) ? 2.0*sqrt(2.0*ip) : 0.5;
    gamm  = pow(qtot, 1.0/3.0)/(2.0*alpha);
  }else{
    alpha = 0.5;
    gamm  = 1.0;
  }

  for(is=0; is<p->nspin; is++){
    gdm = sqrt(sigma[is==0 ? 0 : 2]);

    if(rho[is]>p->threshold && gdm>p->threshold){
      double f;
      
      x =  gdm/pow(rho[is], 4.0/3.0);
      f = -beta*pow(rho[is], 1.0/3.0)*
	x*x/(1.0 + 3.0*beta*x*asinh(gamm*x));
      dedd[is] += f;

    }else if(r > 0.0){
      x = r + (3.0/alpha)*log(2.0*gamm*alpha*pow(qtot, -1.0/3.0));
      /* x = x + pow(qtot*exp(-alpha*r), 1.0/3.0)/(beta*alpha*alpha); */
      dedd[is] -= 1.0/x;
    }
  }
}

func_type func_gga_lb = {
  XC_GGA_XC_LB,
  XC_EXCHANGE_CORRELATION,
  "van Leeuwen & Baerends",
  XC_FAMILY_GGA,
  "R. van Leeuwen and E. J. Baerends, Phys. Rev. A. 49, 2421 (1994)",
  XC_PROVIDES_VXC,
  NULL,
  gga_lb_end,
  NULL,
  NULL /* we can not call this directly */
};

void gga_lb_init(gga_type *p, int nspin, int modified, double threshold)
{
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);

  p->nspin = nspin;
  p->func = &func_gga_lb;
  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_x_init(p->lda_aux, nspin, 3, XC_NON_RELATIVISTIC);
  
  p->modified  = modified;
  p->threshold = threshold;
}
