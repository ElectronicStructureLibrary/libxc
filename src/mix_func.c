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
XC(mix_func)(const XC(func_type) *func,
	     int np, const FLOAT *rho, const FLOAT *sigma,
	     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const XC(func_type) *aux;
  FLOAT *zk_, *vrho_, *vsigma_, *v2rho2_, *v2rhosigma_, *v2sigma2_;
  int ip, ii;

  /* prepare buffers that will hold the results from the individual functionals */
  zk_ = NULL;
  vrho_ = vsigma_ = NULL;
  v2rho2_ = v2rhosigma_ = v2sigma2_ = NULL;

  if(zk != NULL)
    zk_ = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_zk);

  if(vrho != NULL){
    vrho_ = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_vrho);
    if(func->info->family >  XC_FAMILY_LDA){
      vsigma_ = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_vsigma);
    }
  }

  if(v2rho2 != NULL){
    v2rho2_ = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_v2rho2);
    if(func->info->family >  XC_FAMILY_LDA){
      v2rhosigma_ = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_v2rhosigma);
      v2sigma2_   = (FLOAT *) malloc(sizeof(FLOAT)*np*func->n_v2sigma2);
    }
  }

  /* we now add the different components */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    switch(aux->info->family){
    case XC_FAMILY_LDA:
      XC(lda)(aux, np, rho, zk_, vrho_, v2rho2_, NULL);
      break;
    case XC_FAMILY_GGA:
      XC(gga)(aux, np, rho, sigma, zk_, vrho_, vsigma_, v2rho2_, v2rhosigma_, v2sigma2_);
      break;
    }

    if(zk != NULL)
      for(ip = 0; ip < np*func->n_zk; ip++)
	zk[ip] += func->mix_coef[ii] * zk_[ip];

    if(vrho != NULL){
      for(ip = 0; ip < np*func->n_vrho; ip++)
	vrho[ip] += func->mix_coef[ii] * vrho_[ip];

      if(func->info->family > XC_FAMILY_LDA && aux->info->family > XC_FAMILY_LDA)
	for(ip = 0; ip < np*func->n_vsigma; ip++)
	  vsigma[ip] += func->mix_coef[ii] * vsigma_[ip];
    }

    if(v2rho2 != NULL){
      for(ip = 0; ip < np*func->n_v2rho2; ip++)
	v2rho2[ip] += func->mix_coef[ii] * v2rho2_[ip];

      if(func->info->family > XC_FAMILY_LDA && aux->info->family > XC_FAMILY_LDA){
	for(ip = 0; ip < np*func->n_v2rhosigma; ip++)
	  v2rhosigma[ip] += func->mix_coef[ii] * v2rhosigma_[ip];

	for(ip = 0; ip < np*func->n_v2sigma2; ip++)
	  v2sigma2[ip] += func->mix_coef[ii] * v2sigma2_[ip];
      }
    }
  }

  /* deallocate internal buffers */
  if(zk_         != NULL) free(zk_);
  if(vrho_       != NULL) free(vrho_);
  if(vsigma_     != NULL) free(vsigma_);
  if(v2rho2_     != NULL) free(v2rho2_);
  if(v2rhosigma_ != NULL) free(v2rhosigma_);
  if(v2sigma2_   != NULL) free(v2sigma2_);
}
