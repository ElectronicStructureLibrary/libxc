#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void gga_init(gga_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_GGA_X_PBE    ||
	 functional == XC_GGA_C_PBE    ||
	 functional == XC_GGA_XC_LB    ||
	 functional == XC_GGA_X_B88);
  
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals that need it */
  switch(functional){
  case(XC_GGA_X_PBE):
    gga_x_pbe_init(p);
    break;

  case(XC_GGA_C_PBE):
    gga_c_pbe_init(p);
    break;

  case(XC_GGA_XC_LB):
    gga_lb_init(p, nspin, 0, 0.0);
    break;

  case(XC_GGA_X_B88):
    gga_x_b88_init(p);
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
  case(XC_GGA_XC_LB) :
    gga_lb_end(p);
    break;
  case(XC_GGA_X_B88) :
    gga_x_b88_end(p);
    break;
  }
}

/* Some useful formulas:

   sigma_st  = grad rho_s . grad rho_t
   vrho_s    = d e / d rho_s
   vsigma_st = d n*e / d sigma_st

if nspin == 2
   rho   = (rho_u, rho_d)
   sigma = (sigma_uu, sigma_du, sigma_dd)
*/
void gga(gga_type *p, double *rho, double *sigma,
	 double *e, double *vrho, double *vsigma)
{
  double dens;

  assert(p!=NULL);
  
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  
  if(dens <= MIN_DENS){
    int i;

    *e = 0.0;
    for(i=0; i<p->nspin; i++) vrho [i] = 0.0;

    vsigma[0] = 0.0;
    if(p->nspin == XC_POLARIZED){
      vsigma[1] = 0.0; vsigma[2] = 0.0;
    }
    return;
  }
  
  switch(p->func->number){
  case(XC_GGA_X_PBE):
    gga_x_pbe(p, rho, sigma, e, vrho, vsigma);
    break;
    
  case(XC_GGA_C_PBE):
    gga_c_pbe(p, rho, sigma, e, vrho, vsigma);
    break;

  case(XC_GGA_XC_LB): {
    int i;
    *e = 0.0;
    for(i=0; i<3*p->nspin; i++) vsigma[i] = 0.0;
    gga_lb(p, rho, sigma, 0.0, 0.0, 0.0, vrho);
    break;
  }
  case(XC_GGA_X_B88):
    gga_x_b88(p, rho, sigma, e, vrho, vsigma);
  }

}
