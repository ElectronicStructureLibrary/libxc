#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void lca_init(lca_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_LCA_LCH ||
	 functional == XC_LCA_OMC);

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals */
  switch(functional){
  case XC_LCA_LCH:
    lca_lch_init(p);
    break;

  case XC_LCA_OMC:
    lca_omc_init(p);
    break;
  }

}


void lca(lca_type *p, double *rho, double *v, double *e, double *dedd, double *dedv)
{
  int i;
  double dens;

  assert(p!=NULL);

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  
  if(dens <= MIN_DENS){
    *e = 0.0;
    for(i=0; i<  p->nspin; i++) dedd [i] = 0.0;
    for(i=0; i<3*p->nspin; i++) dedv [i] = 0.0;
    return;
  }

  int j;
  double rs, drsdd, vs, vs2, kf, dkfdrs, f, s, dsdrs;
  for(i=0; i<  p->nspin; i++){
    
    rs = RS(rho[i]);
    drsdd = -rs/(3.0*rho[i]);
    vs = sqrt(v _(i, 0)*v _(i, 0) + 
	      v _(i, 1)*v _(i, 1) +
	      v _(i, 2)*v _(i, 2));
    vs2 = vs*vs;
    dkfdrs = pow(9.0*M_PI/4.0, 1.0/3.0);
    kf = dkfdrs*rs;
    f = 24.0*M_PI*M_PI;

    switch(p->func->number){
    case XC_LCA_LCH:
      lca_s_lch(rs, s, dsdrs);
      break;

    case XC_LCA_OMC:
      lca_s_omc(rs, s, dsdrs);
      break;
    }

    *e = kf/f*(s - 1.0)*vs2;
    dedd [i] = drsdd/f*( dkfdrs*(s - 1.0) + kf*dsdrs )*vs2;
    for(j=0; j<3; j++) dedv _(i, j) = 2.0*kf/f*(s - 1.0)*v _(i, j);
  }
  
}
