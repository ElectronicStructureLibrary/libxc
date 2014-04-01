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
#include "funcs_hyb_gga.c"

/* initialization */
int XC(gga_init)(XC(func_type) *func, const XC(func_info_type) *info, int nspin)
{
  assert(func != NULL);

  /* initialize structure */
  func->info   = info;
  func->nspin  = nspin;
  func->params = NULL;
  func->func   = 0;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix_coef   = NULL;
  func->cam_omega = func->cam_alpha = func->cam_beta = 0.0;

  /* initialize spin counters */
  func->n_zk  = 1;
  func->n_rho = func->n_vrho = func->nspin;
  if(func->nspin == XC_UNPOLARIZED){
    func->n_sigma  = func->n_vsigma = 1;
    func->n_v2rho2 = func->n_v2rhosigma  = func->n_v2sigma2 = 1;
    func->n_v3rho3 = func->n_v3rho2sigma = func->n_v3rhosigma2 = func->n_v3sigma3 = 1;
  }else{
    func->n_sigma      = func->n_vsigma = 3;
    func->n_v2rho2     = 3;
    func->n_v2rhosigma = func->n_v2sigma2 = 6;

    func->n_v3rho3      = 4;
    func->n_v3rho2sigma = 9;
    func->n_v3rhosigma2 = 12;
    func->n_v3sigma3    = 10;
  }

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);
  return 0;
}


/* Termination */
void XC(gga_end)(XC(func_type) *func)
{
  assert(func != NULL);

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
    func->n_func_aux = 0;
  }

  if(func->mix_coef != NULL){
    free(func->mix_coef);
    func->mix_coef = NULL;
  }

  /* deallocate any used parameter */
  if(func->params != NULL){
    free(func->params);
    func->params = NULL;
  }
}

/* Some useful formulas:

   sigma_st          = grad rho_s . grad rho_t
   zk                = energy density per unit particle

   vrho_s            = d zk / d rho_s
   vsigma_st         = d n*zk / d sigma_st
   
   v2rho2_st         = d^2 n*zk / d rho_s d rho_t
   v2rhosigma_svx    = d^2 n*zk / d rho_s d sigma_tv
   v2sigma2_stvx     = d^2 n*zk / d sigma_st d sigma_vx

   v3rho3_stv        = d^3 n*zk / d rho_s d rho_t d rho_v
   v3rho2sigma_stvx  = d^3 n*zk / d rho_s d rho_t d sigma_vx
   v3rhosigma2_svxyz = d^3 n*zk / d rho_s d sigma_vx d sigma_yz
   v3sigma3_stvxyz   = d^3 n*zk / d sigma_st d sigma_vx d sigma_yz

if nspin == 2
   rho(2)          = (u, d)
   sigma(3)        = (uu, ud, dd)

   vrho(2)         = (u, d)
   vsigma(3)       = (uu, ud, dd)

   v2rho2(3)       = (u_u, u_d, d_d)
   v2rhosigma(6)   = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
   v2sigma2(6)     = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)

   v3rho3(4)       = (u_u_u, u_u_d, u_d_d, d_d_d)
   v3rho2sigma(9)  = (u_u_uu, u_u_ud, u_u_dd, u_d_uu, u_d_ud, u_d_dd, d_d_uu, d_d_ud, d_d_dd)
   v3rhosigma2(12) = (u_uu_uu, u_uu_ud, u_uu_dd, u_ud_ud, u_ud_dd, u_dd_dd, d_uu_uu, d_uu_ud, d_uu_dd, d_ud_ud, d_ud_dd, d_dd_dd)
   v3sigma(10)     = (uu_uu_uu, uu_uu_ud, uu_uu_dd, uu_ud_ud, uu_ud_dd, uu_dd_dd, ud_ud_ud, ud_ud_dd, ud_dd_dd, dd_dd_dd)
   
*/
void XC(gga)(const XC(func_type) *func, int np, const FLOAT *rho, const FLOAT *sigma,
	     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
	     FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
{
  assert(func != NULL);
  
  /* sanity check */
  if(zk != NULL && !(func->info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc",
	    func->info->name);
    exit(1);
  }

  if(vrho != NULL && !(func->info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc",
	    func->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(func->info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc",
	    func->info->name);
    exit(1);
  }

  if(v3rho3 != NULL && !(func->info->flags & XC_FLAGS_HAVE_KXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of kxc",
	    func->info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL)
    memset(zk, 0, func->n_zk*np*sizeof(FLOAT));

  if(vrho != NULL){
    assert(vsigma != NULL);
    
    memset(vrho,   0, func->n_vrho  *np*sizeof(FLOAT));
    memset(vsigma, 0, func->n_vsigma*np*sizeof(FLOAT));
  }

  if(v2rho2 != NULL){
    assert(v2rhosigma!=NULL && v2sigma2!=NULL);

    memset(v2rho2,     0, func->n_v2rho2    *np*sizeof(FLOAT));
    memset(v2rhosigma, 0, func->n_v2rhosigma*np*sizeof(FLOAT));
    memset(v2sigma2,   0, func->n_v2sigma2  *np*sizeof(FLOAT));
  }

  if(v3rho3 != NULL){
    assert(v3rho2sigma!=NULL && v3rhosigma2!=NULL && v3sigma3!=NULL);

    memset(v3rho3,      0, func->n_v3rho3     *np*sizeof(FLOAT));
    memset(v3rho2sigma, 0, func->n_v3rho2sigma*np*sizeof(FLOAT));
    memset(v3rhosigma2, 0, func->n_v3rhosigma2*np*sizeof(FLOAT));
    memset(v3sigma3,    0, func->n_v3sigma3   *np*sizeof(FLOAT));    
  }

  /* call functional */
  if(func->info->gga != NULL)
    func->info->gga(func, np, rho, sigma, zk, vrho, vsigma, 
		    v2rho2, v2rhosigma, v2sigma2,
		    v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

  if(func->mix_coef != NULL)
    XC(mix_func)(func, np, rho, sigma, NULL, NULL, zk, vrho, vsigma, NULL, NULL,
		 v2rho2, v2sigma2, NULL, NULL, v2rhosigma, NULL, NULL, NULL, NULL, NULL);

}

/* specializations */
/* returns only energy */
inline void 
XC(gga_exc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
	    FLOAT *zk)
{
  XC(gga)(p, np, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

/* returns only potential */
inline void 
XC(gga_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	    FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)(p, np, rho, sigma, NULL, vrho, vsigma, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

/* returns both energy and potential (the most common call usually) */
inline void 
XC(gga_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)(p, np, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

/* returns second derivatives */
inline void 
XC(gga_fxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)(p, np, rho, sigma, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL);
}

/* returns third derivatives */
inline void 
XC(gga_kxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
	    FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
{
  XC(gga)(p, np, rho, sigma, NULL, NULL, NULL, NULL, NULL, NULL, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);
}
