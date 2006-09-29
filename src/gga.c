#include <stdlib.h>
#include <assert.h>

#include "util.h"

extern func_type /* these are the GGA functionals that I know */
  func_gga_x_pbe,
  func_gga_x_b86,
  func_gga_x_b86_r,
  func_gga_x_b88,
  func_gga_c_pbe,
  func_gga_lb;

static func_type *gga_known_funct[] = {
  &func_gga_x_pbe,
  &func_gga_x_b86,
  &func_gga_x_b86_r,
  &func_gga_x_b88,
  &func_gga_c_pbe,
  &func_gga_lb,
  NULL
};


/* initialization */
int gga_init(gga_type *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; gga_known_funct[i]!=NULL; i++){
    if(gga_known_funct[i]->number == functional) break;
  }
  assert(gga_known_funct[i] != NULL);
  if(gga_known_funct[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->func = gga_known_funct[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* see if we need to initialize the functional */
  if(p->func->init != NULL)
    p->func->init(p);
  return 0;
}


void gga_end(gga_type *p)
{
  assert(p != NULL);

  if(p->func->end != NULL)
    p->func->end(p);
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

  assert(p->func!=NULL && p->func->gga!=NULL);
  p->func->gga(p, rho, sigma, e, vrho, vsigma);
}
