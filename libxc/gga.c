#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void gga_init(gga_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_GGA_X_PBE    ||
	 functional == XC_GGA_C_PBE);
  
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals that need it */
  switch(functional){
  case(XC_GGA_X_PBE) :
    gga_x_pbe_init(p);
    break;

  case(XC_GGA_C_PBE) :
    gga_c_pbe_init(p);
    break;
  }
}


void gga_end(gga_type *p)
{
  switch(p->func->number){
  case(XC_GGA_X_PBE) :
  case(XC_GGA_C_PBE) :
    gga_pbe_end(p);
    break;
  }
}


void gga(gga_type *p, double *rho, double *grho,
	 double *e, double *dedd, double *dedgd)
{
  double dens;
  
  assert(p!=NULL);
  
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  
  if(dens <= MIN_DENS){
    int i;
    *e = 0.0;
    for(i=0; i<p->nspin; i++){
      dedd [i] = 0.0;
      dedgd[i] = 0.0;
    }
    return;
  }
  
  switch(p->func->number){
  case(XC_GGA_X_PBE):
    gga_x_pbe(p, rho, grho, e, dedd, dedgd);
    break;

  case(XC_GGA_C_PBE):
    gga_c_pbe(p, rho, grho, e, dedd, dedgd);
    break;
  }

}
