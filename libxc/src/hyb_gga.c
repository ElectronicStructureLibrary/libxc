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
void XC(hyb_gga)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  FLOAT zk_, vrho_[2], vsigma_[3], v2rho2_[6], v2rhosigma_[6], v2sigma2_[6];
  FLOAT *pzk, *pvrho, *pv2rho2;
  int ii, is, js;

  assert(p!=NULL && p->info!=NULL);
  
  if(!XC(gga_input_init)(p->info, p->nspin, rho, zk, vrho, vsigma,
			 v2rho2, v2rhosigma, v2sigma2)) return; 

  /* hybrid may want to add some term */
  if(p->info->gga!=NULL)
    p->info->gga(p, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

  pzk     = (zk     == NULL) ? NULL : &zk_;
  pvrho   = (vrho   == NULL) ? NULL : vrho_;
  pv2rho2 = (v2rho2 == NULL) ? NULL : v2rho2_;

  /* we now add the LDA components */
  for(ii=0; ii<p->lda_n; ii++){
    XC(lda)(p->lda_aux[ii], rho, pzk, pvrho, pv2rho2, NULL);

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

  /* and the GGA components */
  for(ii=0; ii<p->gga_n; ii++){
    XC(gga)(p->gga_aux[ii], rho, sigma, pzk, pvrho, vsigma_, pv2rho2,  v2rhosigma_, v2sigma2_);

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

}

/* especializations */
inline void 
XC(hyb_gga_exc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		FLOAT *zk)
{
  XC(hyb_gga)(p, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(hyb_gga_vxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(hyb_gga)(p, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL);
}

inline void 
XC(hyb_gga_fxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(hyb_gga)(p, rho, sigma, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2);
}


/*****************************************************/
FLOAT XC(hyb_gga_exx_coef)(XC(hyb_gga_type) *p)
{
  assert(p!=NULL);

  return p->exx_coef;
}
