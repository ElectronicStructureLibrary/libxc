#include <stdlib.h>
#include <assert.h>

#include "util.h"

extern xc_func_info_type /* these are the GGA functionals that I know */
  func_info_gga_x_pbe,
  func_info_gga_x_pbe_r,
  func_info_gga_x_b86,
  func_info_gga_x_b86_r,
  func_info_gga_x_b86_mgc,
  func_info_gga_x_b88,
  func_info_gga_x_g96,
  func_info_gga_c_pbe,
  func_info_gga_c_lyp,
  func_info_gga_lb;

xc_func_info_type *gga_known_funct[] = {
  &func_info_gga_x_pbe,
  &func_info_gga_x_pbe_r,
  &func_info_gga_x_b86,
  &func_info_gga_x_b86_r,
  &func_info_gga_x_b86_mgc,
  &func_info_gga_x_b88,
  &func_info_gga_x_g96,
  &func_info_gga_c_pbe,
  &func_info_gga_c_lyp,
  &func_info_gga_lb,
  NULL
};


/* initialization */
int xc_gga_init(xc_gga_type *p, int functional, int nspin)
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
  p->info = gga_known_funct[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}


/* Termination */
void xc_gga_end(xc_gga_type *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);
}


/* Some useful formulas:

   sigma_st  = grad rho_s . grad rho_t
   vrho_s    = d e / d rho_s
   vsigma_st = d n*e / d sigma_st

if nspin == 2
   rho   = (rho_u, rho_d)
   sigma = (sigma_uu, sigma_du, sigma_dd)
*/
void xc_gga(xc_gga_type *p, double *rho, double *sigma,
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

  assert(p->info!=NULL && p->info->gga!=NULL);
  p->info->gga(p, rho, sigma, e, vrho, vsigma);
}
