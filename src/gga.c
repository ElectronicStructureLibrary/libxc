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
#include <string.h>
#include <assert.h>

#include "util.h"
#include "funcs_gga.c"


/* initialization */
int XC(gga_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin)
{
  XC(gga_type) *func;

  assert(p != NULL && p->gga != NULL);
  func = p->gga;

  /* initialize structure */
  func->info   = info;
  func->nspin  = nspin;
  func->params = NULL;
  func->func   = 0;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix        = NULL;

  /* initialize spin counters */
  func->n_zk  = 1;
  func->n_rho = func->n_vrho = func->nspin;
  if(nspin == XC_UNPOLARIZED){
    func->n_sigma  = func->n_vsigma = 1;
    func->n_v2rho2 = func->n_v2rhosigma = func->n_v2sigma2 =1;
  }else{
    func->n_sigma      = func->n_vsigma = func->n_v2rho2 = 3;
    func->n_v2rhosigma = func->n_v2sigma2 = 6;
  }

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);
  return 0;
}


/* Termination */
void XC(gga_end)(XC(func_type) *p)
{
  XC(gga_type) *func;

  assert(p != NULL && p->gga != NULL);
  func = p->gga;

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(func);

  /* terminate any auxiliary functional */
  if(func->n_func_aux > 0){
    int ii;

    for(ii=0; ii<func->n_func_aux; ii++){
      XC(func_end)(func->func_aux[ii]);
      free(func->func_aux[ii]);
    }
    free(func->func_aux);
  }

  /* deallocate any used parameter */
  if(func->params != NULL){
    free(func->params);
    func->params = NULL;
  }

  if(func->mix != NULL){
    XC(mix_func_free)(func->mix);
    free(func->mix);
    func->mix = NULL;
  }
}

/* Some useful formulas:

   sigma_st       = grad rho_s . grad rho_t
   zk             = energy density per unit particle

   vrho_s         = d zk / d rho_s
   vsigma_st      = d n*zk / d sigma_st
   
   v2rho2_st      = d^2 n*zk / d rho_s d rho_t
   v2rhosigma_svx = d^2 n*zk / d rho_s d sigma_tv
   v2sigma2_stvx  = d^2 n*zk / d sigma_st d sigma_vx

if nspin == 2
   rho(2)        = (u, d)
   sigma(3)      = (uu, du, dd)

   vrho(2)       = (u, d)
   vsigma(3)     = (uu, du, dd)

   v2rho2(3)     = (uu, du, dd)
   v2rhosigma(6) = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
   v2sigma2(6)   = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)
*/
void XC(gga)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga_type) *func;

  assert(p != NULL && p->gga != NULL);
  func = p->gga;
  
  /* sanity check */
  if(zk != NULL && !(func->info->provides & XC_PROVIDES_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc",
	    func->info->name);
    exit(1);
  }

  if(vrho != NULL && !(func->info->provides & XC_PROVIDES_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc",
	    func->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(func->info->provides & XC_PROVIDES_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc",
	    func->info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL)
    memset(zk, 0.0, func->n_zk*np*sizeof(FLOAT));

  if(vrho != NULL){
    assert(vsigma != NULL);
    
    memset(vrho,   0.0, func->n_vrho  *np*sizeof(FLOAT));
    memset(vsigma, 0.0, func->n_vsigma*np*sizeof(FLOAT));
  }

  if(v2rho2 != NULL){
    assert(v2rhosigma!=NULL && v2sigma2!=NULL);

    memset(v2rho2,     0.0, func->n_v2rho2    *np*sizeof(FLOAT));
    memset(v2rhosigma, 0.0, func->n_v2rhosigma*np*sizeof(FLOAT));
    memset(v2sigma2,   0.0, func->n_v2sigma2  *np*sizeof(FLOAT));
  }

  /* call functional */
  if(func->info->gga != NULL)
    func->info->gga(func, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

  if(func->mix != NULL)
    XC(mix_func)(func->mix, np, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

/* especializations */
/* returns only energy */
inline void 
XC(gga_exc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
	    FLOAT *zk)
{
  XC(gga)(p, np, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL);
}

/* returns only potential */
inline void 
XC(gga_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	    FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)(p, np, rho, sigma, NULL, vrho, vsigma, NULL, NULL, NULL);
}

/* returns both energy and potential (the most common call usually) */
inline void 
XC(gga_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)(p, np, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL);
}

/* returns second derivatives */
inline void 
XC(gga_fxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)(p, np, rho, sigma, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2);
}
