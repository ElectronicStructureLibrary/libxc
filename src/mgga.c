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
#include "funcs_mgga.c"

/* initialization */
int XC(mgga_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin)
{
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

  /* initialize structure */
  func->info   = info;
  func->nspin  = nspin;
  func->params = NULL;
  func->func   = 0;

  func->n_func_aux = 0;
  func->func_aux   = NULL;

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);
  return 0;
}


void XC(mgga_end)(XC(func_type) *p)
{
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(p);

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
}


void 
XC(mgga)(const XC(func_type) *p, 
	 const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  FLOAT dens;
  int i, n;
  XC(mgga_type) *func;

  assert(p != NULL && p->mgga != NULL);
  func = p->mgga;

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
  if(zk != NULL){
    *zk = 0.0;
  }

  if(vrho != NULL){
    assert(vsigma != NULL);

    for(i=0; i<func->nspin; i++){
      vrho[i] = 0.0;
      vtau[i] = 0.0;
    }

    n = (func->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++)
      vsigma[i] = 0.0;
  }

  if(v2rho2 != NULL){
    assert(v2rhosigma!=NULL && v2sigma2!=NULL && v2rhotau!=NULL && v2tausigma!=NULL && v2tau2!=NULL);

    n = (func->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++){
      v2rho2[i] = 0.0;
      v2tau2[i] = 0.0;
    }

    n = (func->nspin == XC_UNPOLARIZED) ? 1 : 4;
    for(i=0; i<n; i++){
      v2rhotau[i] = 0.0;
    }

    n = (func->nspin == XC_UNPOLARIZED) ? 1 : 6;
    for(i=0; i<n; i++){
      v2rhosigma[i] = 0.0;
      v2tausigma[i] = 0.0;
    }
  }

  /* check if density is larger than threshold */
  dens = rho[0];
  if(func->nspin == XC_POLARIZED) dens += rho[1];
  if(dens <= MIN_DENS) return;

  /* call functional */
  assert(func->info->mgga != NULL);
  func->info->mgga(func, rho, sigma, lapl_rho, tau, zk, vrho, vsigma, vlapl_rho, vtau, 
		v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}

/* especializations */
inline void 
XC(mgga_exc)(const XC(func_type) *p, const FLOAT *rho, 
	     const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *zk)
{
  XC(mgga)(p, rho, sigma, tau, zk, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_exc_vxc)(const XC(func_type) *p, const FLOAT *rho,
		 const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga)(p, rho, sigma, lapl_rho, tau, zk, vrho, vsigma, vlapl_rho, vtau, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_vxc)(const XC(func_type) *p, const FLOAT *rho, 
	     const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga)(p, rho, sigma, lapl_rho, tau, NULL, vrho, vsigma, vlapl_rho, vtau, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(mgga_fxc)(const XC(func_type) *p, const FLOAT *rho, const FLOAT *sigma,
	     const FLOAT *lapl_rho, const FLOAT *tau,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga)(p, rho, sigma, lapl_rho, tau, NULL, NULL, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}


