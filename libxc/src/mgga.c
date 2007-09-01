#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void xc_mgga_init(xc_mgga_type *p, int functional, int nspin)
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
  case(XC_MGGA_C_TPSS):
    mgga_c_tpss_init(p);
    break;
  }
}


void xc_mgga_end(xc_mgga_type *p)
{
  switch(p->info->number){
  case(XC_MGGA_X_TPSS) :
    mgga_x_tpss_end(p);
    break;
  case(XC_MGGA_C_TPSS) :
    mgga_c_tpss_end(p);
    break;
  }
}


void xc_mgga(xc_mgga_type *p, double *rho, double *grho, double *tau,
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
  
  switch(p->info->number){
  case(XC_MGGA_X_TPSS):
    mgga_x_tpss(p, rho, grho, tau, e, dedd, dedgd, dedtau);
    break;

  case(XC_MGGA_C_TPSS):
    mgga_c_tpss(p, rho, grho, tau, e, dedd, dedgd, dedtau);
    break;
  }

}

void xc_mgga_sp(xc_mgga_type *p, float *rho, float *grho, float *tau,
		float *e, float *dedd, float *dedgd, float *dedtau){

  double drho[2], dgrho[6], dtau[2];
  double de[1], ddedd[2], ddedgd[6], ddedtau[2];
  int ii;

  for(ii=0; ii < p->nspin; ii++) drho[ii] = rho[ii];
  for(ii=0; ii < 3*p->nspin; ii++) dgrho[ii] = grho[ii];
  for(ii=0; ii < p->nspin; ii++) dtau[ii] = tau[ii];

  xc_mgga(p, drho, dgrho, dtau,
	  de, ddedd, ddedgd, ddedtau);
  
  e[0] = de[0];
  for(ii=0; ii < p->nspin; ii++) dedd[ii] = ddedd[ii];
  for(ii=0; ii < 3*p->nspin; ii++) dedgd[ii] = ddedgd[ii];
  for(ii=0; ii < p->nspin; ii++) dedtau[ii] = ddedtau[ii];

}

