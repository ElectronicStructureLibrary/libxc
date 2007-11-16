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

/* Note: Do not forget to add a correlation (LDA) functional to the LB94 */
#define XC_GGA_XC_LB 160 /* van Leeuwen & Baerends */

typedef struct{
  int    modified;  /* shall we used a modified version */
  double threshold; /* when to start using the analitic form */
  double ip;        /* ionization potential of the species */
  double qtot;      /* total charge in the region */

  double alpha;     /* the parameters of LB94 */
  double gamm;
} gga_xc_lb_params;

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

void gga_lb_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  assert(p->params == NULL);

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_x_init(p->lda_aux, p->nspin, 3, XC_NON_RELATIVISTIC);

  p->params = malloc(sizeof(gga_xc_lb_params));
  xc_gga_lb_set_params(p, 0, 0.0, 0.0, 0.0);
}


void gga_lb_end(void *p_)
{
  xc_gga_type *p = p_;

  free(p->lda_aux);
  free(p->params);
}


void xc_gga_lb_set_params(xc_gga_type *p, int modified, double threshold, double ip, double qtot)
{
  gga_xc_lb_params *params;

  fflush(stdout);
  assert(p->params != NULL);
  params = (gga_xc_lb_params *) (p->params);

  params->modified  = modified;
  params->threshold = threshold;
  params->ip        = ip;
  params->qtot      = qtot;

  if(params->modified){
    params->alpha = (params->ip > 0.0) ? 2.0*sqrt(2.0*params->ip) : 0.5;
    params->gamm  = pow(params->qtot, 1.0/3.0)/(2.0*params->alpha);
  }else{
    params->alpha = 0.5;
    params->gamm  = 1.0;
  }
}


void xc_gga_lb_set_params_sp(xc_gga_type *p, int modified, float threshold, float ip, float qtot)
{
  xc_gga_lb_set_params(p, modified, (double)threshold, (double)ip, (double)qtot);
}


void xc_gga_lb_modified(xc_gga_type *p, double *rho, double *sigma, double r, double *dedd)
{
  int is;
  double gdm, x;
  gga_xc_lb_params *params = (gga_xc_lb_params *) (p->params);

  static const double beta = 0.05;

  xc_lda_vxc(p->lda_aux, rho, &x, dedd);

  for(is=0; is<p->nspin; is++){
    gdm = sqrt(sigma[is==0 ? 0 : 2]);

    if(params->modified == 0 || 
       (rho[is] > params->threshold && gdm > params->threshold)){
      double f;
      
      if(rho[is] <= MIN_DENS) continue;

      x =  gdm/pow(rho[is], 4.0/3.0);

      /* dirty fix, destroys -1/r, but results will not run wild */
      if(params->modified == 0 && x > 500.0) continue;

      /* the actual functional */
      f = -beta*pow(rho[is], 1.0/3.0)*
	x*x/(1.0 + 3.0*beta*x*asinh(params->gamm*x));
      dedd[is] += f;

    }else if(r > 0.0){
      /* the aymptotic expansion of LB94 */
      x = r + (3.0/params->alpha)*
	log(2.0*params->gamm * params->alpha * pow(params->qtot, -1.0/3.0));

      /* x = x + pow(qtot*exp(-alpha*r), 1.0/3.0)/(beta*alpha*alpha); */

      dedd[is] -= 1.0/x;
    }
  }
}


void xc_gga_lb_modified_sp(xc_gga_type *p, float *rho, float *sigma, float r, float *dedd){

  double drho[2], dsigma[6];
  int ii, nsig;
  double ddedd[2];

  for(ii=0; ii < p->nspin; ii++) drho[ii] = rho[ii];
  
  nsig = (p->nspin == XC_POLARIZED) ? 1 : 3;
  for(ii=0; ii < nsig; ii++) dsigma[ii] = sigma[ii];
  
  xc_gga_lb_modified(p, drho, dsigma, (double) r, ddedd);
  
  for(ii=0; ii < p->nspin; ii++) dedd[ii] = (float)ddedd[ii];
}


static void gga_xc_lb(void *p_, double *rho, double *sigma,
                      double *e, double *vrho, double *vsigma)
{
  xc_gga_lb_modified((xc_gga_type *)p_, rho, sigma, 0.0, vrho);
}


xc_func_info_type func_info_gga_xc_lb = {
  XC_GGA_XC_LB,
  XC_EXCHANGE_CORRELATION,
  "van Leeuwen & Baerends",
  XC_FAMILY_GGA,
  "R van Leeuwen and EJ Baerends, Phys. Rev. A. 49, 2421 (1994)",
  XC_PROVIDES_VXC,
  gga_lb_init,
  gga_lb_end,
  NULL,
  gga_xc_lb
};

