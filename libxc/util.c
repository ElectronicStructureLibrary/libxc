#include <assert.h>

#include "util.h"

/* this function converts the spin-density into total density and
	 relative magnetization */
void rho2dzeta(int nspin, double *rho, double *d, double *zeta)
{
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  
  if(nspin==XC_UNPOLARIZED){
    *d    = max(MIN_DENS, rho[0]);
    *zeta = 0.0;
  }else{
    *d    = max(MIN_DENS, rho[0]+rho[1]);
    *zeta = (*d > MIN_DENS) ? (rho[0]-rho[1])/(*d) : 0.0;
  }
}
