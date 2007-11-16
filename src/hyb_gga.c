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

#include <stdlib.h>
#include <assert.h>

#include "util.h"
#include "funcs_hyb_gga.c"

/* initialization */
/*****************************************************/
int xc_hyb_gga_init(xc_hyb_gga_type *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; hyb_gga_known_funct[i]!=NULL; i++){
    if(hyb_gga_known_funct[i]->number == functional) break;
  }
  if(hyb_gga_known_funct[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->params = NULL;
  p->info = hyb_gga_known_funct[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  p->lda_n = 0;
  p->gga_n = 0;
  p->exx_coef = 0.0;

  /* we always need to initialize the functional */
  p->info->init(p);
  return 0;
}


void xc_hyb_gga_alloc(xc_hyb_gga_type *p)
{
  int i;

  if(p->lda_n > 0){
    p->lda_coef = (double *)       malloc(p->lda_n * sizeof(double));
    p->lda_aux  = (xc_lda_type **) malloc(p->lda_n * sizeof(xc_lda_type *));
    for(i=0; i<p->lda_n; i++)
      p->lda_aux[i] = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  }

  if(p->gga_n > 0){
    p->gga_coef = (double *)       malloc(p->gga_n * sizeof(double));
    p->gga_aux  = (xc_gga_type **) malloc(p->gga_n * sizeof(xc_gga_type *));
    for(i=0; i<p->gga_n; i++)
      p->gga_aux[i] = (xc_gga_type *) malloc(sizeof(xc_gga_type));
  }
}

/* Termination */
/*****************************************************/
void xc_hyb_gga_end(xc_hyb_gga_type *p)
{
  int i;

  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);

  /* free the LDA components */
  if(p->lda_n > 0){
    for(i=0; i<p->lda_n; i++)
      free(p->lda_aux[i]);
    free(p->lda_aux);
    free(p->gga_coef);
  }

  /* free the GGA components */
  if(p->gga_n > 0){
    for(i=0; i<p->gga_n; i++){
      xc_gga_end(p->gga_aux[i]);
      free(p->gga_aux[i]);
    }
    free(p->gga_aux);
    free(p->gga_coef);
  }
}


/*****************************************************/
void xc_hyb_gga(xc_hyb_gga_type *p, double *rho, double *sigma,
		double *e, double *vrho, double *vsigma)
{
  double dens;
  int i;

  assert(p!=NULL);
  
  double e1, vrho1[2], vsigma1[3];
  int ii, is, js = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;

  /* initialize all to zero */
  *e = 0.0;
  for(is=0; is<p->nspin; is++) vrho[is]   = 0.0;
  for(is=0; is<js;       is++) vsigma[is] = 0.0;

  /* hybrid may want to add some term */
  if(p->info!=NULL && p->info->gga!=NULL){
    p->info->gga(p, rho, sigma, e, vrho, vsigma);
  }

  /* we now add the LDA components */
  for(ii=0; ii<p->lda_n; ii++){
    xc_lda_vxc(p->lda_aux[ii], rho, &e1, vrho1);

    *e += p->lda_coef[ii] * e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += p->lda_coef[ii] * vrho1[is];
  }

  /* and the GGA components */
  for(ii=0; ii<p->gga_n; ii++){
    xc_gga(p->gga_aux[ii], rho, sigma, &e1, vrho1, vsigma1);

    *e += p->gga_coef[ii] * e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += p->gga_coef[ii] * vrho1[is];

    for(is=0; is<js; is++)
      vsigma[is] += p->gga_coef[ii] * vsigma1[is];
  }

}


void xc_hyb_gga_sp(xc_hyb_gga_type *p, float *rho, float *sigma,
		   float *e, float *vrho, float *vsigma)
{
  double drho[2], dsigma[6];
  double de[1], dvrho[2], dvsigma[6];
  int ii, nsig;

  nsig = (p->nspin == XC_POLARIZED) ? 1 : 3;
  for(ii=0; ii < p->nspin; ii++) drho[ii]   = (double) rho[ii];
  for(ii=0; ii < nsig;     ii++) dsigma[ii] = (double) sigma[ii];

  xc_hyb_gga(p, drho, dsigma, de, dvrho, dvsigma);
  
  e[0] = (float)de[0];
  for(ii=0; ii < p->nspin; ii++) vrho[ii]   = (float)dvrho[ii];
  for(ii=0; ii < nsig    ; ii++) vsigma[ii] = (float)dvsigma[ii];
}


/*****************************************************/
double xc_hyb_gga_exx_coef(xc_hyb_gga_type *p)
{
  assert(p!=NULL);

  return p->exx_coef;
}


float xc_hyb_gga_exx_coef_sp(xc_hyb_gga_type *p)
{
  assert(p!=NULL);

  return (float)p->exx_coef;
}
