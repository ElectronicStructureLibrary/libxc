#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void mgga_init(mgga_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_MGGA_X_TPSS    ||
	 functional == XC_MGGA_C_TPSS);
  
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals that need it */
  switch(functional){
  case(XC_MGGA_X_TPSS):
    mgga_x_tpss_init(p);
    break;
  }
}


void mgga_end(mgga_type *p)
{
  switch(p->func->number){
  case(XC_MGGA_X_TPSS) :
    mgga_x_tpss_end(p);
    break;
  }
}


void mgga(mgga_type *p, double *rho, double *grho, double *tau,
	  double *e, double *dedd, double *dedgd, double *dedtau)

{
  double dens;

  assert(p!=NULL);
  
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  
  if(dens <= MIN_DENS){
    int i;
    *e = 0.0;
    for(i=0; i<  p->nspin; i++){
      dedd  [i] = 0.0;
      dedtau[i] = 0.0;
    }
    for(i=0; i<3*p->nspin; i++) dedgd[i] = 0.0;
    return;
  }
  
  switch(p->func->number){
  case(XC_MGGA_X_TPSS):
    mgga_x_tpss(p, rho, grho, tau, e, dedd, dedgd, dedtau);
    break;
  }

}
