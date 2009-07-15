/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"


/*****************************************************/
void 
XC(mix_func_init)(XC(mix_func_type) *p, int level, int nspin)
{
  p->level = level;
  p->nspin = nspin;

  p->lda_n    = p->gga_n    = p->mgga_n    = 0;
  p->lda_coef = p->gga_coef = p->mgga_coef = NULL;

  p->lda_mix   = NULL;
  p->gga_mix   = NULL;
  p->mgga_mix  = NULL;
}


/*****************************************************/
void 
XC(mix_func_alloc)(XC(mix_func_type) *p)
{
  if(p->lda_n > 0){
    p->lda_coef = (FLOAT *)        malloc(p->lda_n * sizeof(FLOAT));
    p->lda_mix  = (XC(lda_type) *) malloc(p->lda_n * sizeof(XC(lda_type)));
  }

  if(p->gga_n > 0){
    p->gga_coef = (FLOAT *)        malloc(p->gga_n * sizeof(FLOAT));
    p->gga_mix  = (XC(gga_type) *) malloc(p->gga_n * sizeof(XC(gga_type)));
  }
}


/*****************************************************/
void 
XC(mix_func_free)(XC(mix_func_type) *p)
{
  int i;

  if(p->lda_n > 0){
    free(p->lda_mix);
    free(p->lda_coef);
  }

  if(p->gga_n > 0){
    for(i=0; i<p->gga_n; i++)
      XC(gga_end)(&p->gga_mix[i]);
    
    free(p->gga_mix);
    free(p->gga_coef);
  }

  if(p->mgga_n > 0){
    for(i=0; i<p->mgga_n; i++)
      XC(mgga_end)(&p->mgga_mix[i]);
    
    free(p->mgga_mix);
    free(p->mgga_coef);
  }

  XC(mix_func_init)(p, -1, -1);
}


/*****************************************************/
void 
XC(mix_func)(XC(mix_func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  FLOAT zk_, vrho_[2], vsigma_[3], v2rho2_[6], v2rhosigma_[6], v2sigma2_[6];
  FLOAT *pzk, *pvrho, *pv2rho2;
  int ip, ii, is, js;

  exit(1); /* NP HAS TO BE TAKEN INTO ACCOUNT */

  pzk     = (zk     == NULL) ? NULL : &zk_;
  pvrho   = (vrho   == NULL) ? NULL : vrho_;
  pv2rho2 = (v2rho2 == NULL) ? NULL : v2rho2_;

  /* we now add the LDA components */
  for(ip=0; ip<np; ip++){
    for(ii=0; ii<p->lda_n; ii++){
      XC(lda)(&p->lda_mix[ii], 1, rho, pzk, pvrho, pv2rho2, NULL);

      if(zk != NULL)
	*zk += p->lda_coef[ii] * zk_;

      if(vrho != NULL)
	for(is=0; is<p->nspin; is++)
	  vrho[is] += p->lda_coef[ii] * vrho_[is];

      if(v2rho2 != NULL){
	js = (p->nspin == XC_POLARIZED) ? 6 : 1;
	for(is=0; is<js; is++)
	  v2rho2[is] += p->lda_coef[ii] * v2rho2_[is];
      }
    }

    if(p->level == XC_FAMILY_LDA) goto end_ip_loop;
    
    /* and the GGA components */
    for(ii=0; ii<p->gga_n; ii++){
      XC(gga)(&p->gga_mix[ii], 1, rho, sigma, pzk, pvrho, vsigma_, pv2rho2,  v2rhosigma_, v2sigma2_);
      
      if(zk != NULL)
	*zk += p->gga_coef[ii] * zk_; 

      if(vrho != NULL){
	for(is=0; is<p->nspin; is++)
	  vrho[is] += p->gga_coef[ii] * vrho_[is];
	
	js = (p->nspin == XC_POLARIZED) ? 3 : 1;
	for(is=0; is<js; is++)
	  vsigma[is] += p->gga_coef[ii] * vsigma_[is];
      }
      
      if(v2rho2 != NULL){
	js = (p->nspin == XC_POLARIZED) ? 6 : 1;
	for(is=0; is<js; is++){
	  v2rho2[is]     += p->gga_coef[ii] * v2rho2_[is];
	  v2rhosigma[is] += p->gga_coef[ii] * v2rhosigma_[is];
	  v2sigma2[is]   += p->gga_coef[ii] * v2sigma2_[is];
	}
      }
    }

    if(p->level == XC_FAMILY_GGA) goto end_ip_loop;

  end_ip_loop:

    /* increment pointers */
    rho += p->lda_mix[0].n_rho;
    if(zk     != NULL) zk     += p->lda_mix[0].n_zk;
    if(vrho   != NULL) vrho   += p->lda_mix[0].n_vrho;
    if(v2rho2 != NULL) v2rho2 += p->lda_mix[0].n_v2rho2;

    if(p->level == XC_FAMILY_LDA) continue;

    sigma += p->gga_mix[0].n_sigma;
    if(vrho != NULL) vsigma += p->gga_mix[0].n_vsigma;
    if(v2rho2 != NULL){
      v2rhosigma += p->gga_mix[0].n_v2rhosigma;
      v2sigma2   += p->gga_mix[0].n_v2sigma2;
    }
  }
}
