/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_GGA_XC_LB 160 /* van Leeuwen & Baerends */

typedef struct{
  int modified;
  double threshold;
} gga_xc_lb_params;

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

void gga_lb_end(void *p_)
{
  xc_gga_type *p = p_;

  free(p->lda_aux);
  free(p->params);
}

void xc_gga_lb(xc_gga_type *p, double *rho, double *sigma, double r, double ip, double qtot,
	    double *dedd)
{
  int is;
  double alpha, gdm, gamm, x;
  gga_xc_lb_params *params = (gga_xc_lb_params *) (p->params);

  static const double beta = 0.05;

  xc_lda_vxc(p->lda_aux, rho, &x, dedd);
  if(params->modified){
    alpha = (ip > 0.0) ? 2.0*sqrt(2.0*ip) : 0.5;
    gamm  = pow(qtot, 1.0/3.0)/(2.0*alpha);
  }else{
    alpha = 0.5;
    gamm  = 1.0;
  }

  for(is=0; is<p->nspin; is++){
    gdm = sqrt(sigma[is==0 ? 0 : 2]);

    if(rho[is] > params->threshold && gdm > params->threshold){
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

xc_func_info_type func_info_gga_xc_lb = {
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

void xc_gga_lb_init(xc_gga_type *p, int nspin, int modified, double threshold)
{
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);

  p->nspin = nspin;
  p->info = &func_info_gga_xc_lb;
  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_x_init(p->lda_aux, nspin, 3, XC_NON_RELATIVISTIC);
  
  p->params = malloc(sizeof(gga_xc_lb_params));

  {
    gga_xc_lb_params *params = (gga_xc_lb_params *) (p->params);

    params->modified  = modified;
    params->threshold = threshold;
  }
}

void xc_gga_lb_sp(xc_gga_type *p, float *rho, float *sigma, float r, float ip, float qtot,
		  float *dedd){

  double drho[2], dsigma[6];
  int ii;
  double ddedd[2];

  for(ii=0; ii < p->nspin; ii++) drho[ii] = rho[ii];
  for(ii=0; ii < 3*p->nspin; ii++) dsigma[ii] = sigma[ii];

  xc_gga_lb(p, drho, dsigma, (double) r, (double) ip, (double) qtot, ddedd);

  for(ii=0; ii < p->nspin; ii++) dedd[ii] = ddedd[ii];
}
