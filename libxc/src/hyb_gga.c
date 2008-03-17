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
#include "funcs_hyb_gga.c"

/* initialization */
/*****************************************************/
int XC(hyb_gga_init)(XC(hyb_gga_type) *p, int functional, int nspin)
{
  int i;

  assert(p!=NULL);

  /* let us first find out if we know the functional */
  for(i=0; XC(hyb_gga_known_funct)[i]!=NULL; i++){
    if(XC(hyb_gga_known_funct)[i]->number == functional) break;
  }
  if(XC(hyb_gga_known_funct)[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->params = NULL;
  p->info = XC(hyb_gga_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  p->lda_n = 0;
  p->gga_n = 0;
  p->exx_coef = 1.0;

  /* we always need to initialize the functional */
  assert(p->info->init != NULL);
  p->info->init(p);
  return 0;
}


void XC(hyb_gga_alloc)(XC(hyb_gga_type) *p)
{
  int i;

  if(p->lda_n > 0){
    p->lda_coef = (FLOAT *)         malloc(p->lda_n * sizeof(FLOAT));
    p->lda_aux  = (XC(lda_type) **) malloc(p->lda_n * sizeof(XC(lda_type) *));
    for(i=0; i<p->lda_n; i++)
      p->lda_aux[i] = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  }

  if(p->gga_n > 0){
    p->gga_coef = (FLOAT *)         malloc(p->gga_n * sizeof(FLOAT));
    p->gga_aux  = (XC(gga_type) **) malloc(p->gga_n * sizeof(XC(gga_type) *));
    for(i=0; i<p->gga_n; i++)
      p->gga_aux[i] = (XC(gga_type) *) malloc(sizeof(XC(gga_type)));
  }
}


/* Termination */
/*****************************************************/
void XC(hyb_gga_end)(XC(hyb_gga_type) *p)
{
  int i;

  assert(p!=NULL);

  if(p->info->end != NULL)
    p->info->end(p);

  /* free the LDA components */
  if(p->lda_n > 0){
    for(i=0; i<p->lda_n; i++)
      free(p->lda_aux[i]);
    free(p->lda_aux);
    free(p->lda_coef);
  }

  /* free the GGA components */
  if(p->gga_n > 0){
    for(i=0; i<p->gga_n; i++){
      XC(gga_end)(p->gga_aux[i]);
      free(p->gga_aux[i]);
    }
    free(p->gga_aux);
    free(p->gga_coef);
  }
}


/*****************************************************/
void XC(hyb_gga)(XC(hyb_gga_type) *p, FLOAT *rho, FLOAT *sigma,
		FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  FLOAT e1, vrho1[2], vsigma1[3];
  int ii, is, js = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;

  assert(p!=NULL);
  
  /* initialize all to zero */
  *e = 0.0;
  for(is=0; is<p->nspin; is++) vrho[is]   = 0.0;
  for(is=0; is<js;       is++) vsigma[is] = 0.0;

  /* hybrid may want to add some term */
  if(p->info!=NULL && p->info->gga!=NULL){
    p->info->gga(p, rho, sigma, e, vrho, vsigma, NULL, NULL, NULL);
  }

  /* we now add the LDA components */
  for(ii=0; ii<p->lda_n; ii++){
    XC(lda_vxc)(p->lda_aux[ii], rho, &e1, vrho1);

    *e += p->lda_coef[ii] * e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += p->lda_coef[ii] * vrho1[is];
  }

  /* and the GGA components */
  for(ii=0; ii<p->gga_n; ii++){
    XC(gga_vxc)(p->gga_aux[ii], rho, sigma, &e1, vrho1, vsigma1);

    *e += p->gga_coef[ii] * e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += p->gga_coef[ii] * vrho1[is];

    for(is=0; is<js; is++)
      vsigma[is] += p->gga_coef[ii] * vsigma1[is];
  }

}


/*****************************************************/
FLOAT XC(hyb_gga_exx_coef)(XC(hyb_gga_type) *p)
{
  assert(p!=NULL);

  return p->exx_coef;
}

