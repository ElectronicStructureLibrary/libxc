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
#include "funcs_gga.c"


/* initialization */
int XC(gga_init)(XC(gga_type) *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; XC(gga_known_funct)[i]!=NULL; i++){
    if(XC(gga_known_funct)[i]->number == functional) break;
  }
  if(XC(gga_known_funct)[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->mix    = NULL;
  p->params = NULL;
  p->info = XC(gga_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}


/* Termination */
void XC(gga_end)(XC(gga_type) *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);

  if(p->mix != NULL)
    XC(mix_func_free)(p->mix);
  free(p->mix); p->mix = NULL;
}


inline int
XC(gga_input_init)(const XC(func_info_type) *info, int nspin, const FLOAT *rho,
		   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  FLOAT dens;
  int i, n;

  /* sanity check */
  if(zk != NULL && !(info->provides & XC_PROVIDES_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc",
	    info->name);
    exit(1);
  }

  if(vrho != NULL && !(info->provides & XC_PROVIDES_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc",
	    info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(info->provides & XC_PROVIDES_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc",
	    info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL){
    assert(info->provides & XC_PROVIDES_EXC);
    *zk = 0.0;
  }

  if(vrho != NULL){
    assert(info->provides & XC_PROVIDES_VXC);
    assert(vsigma != NULL);
    
    for(i=0; i<nspin; i++) vrho [i] = 0.0;
    
    n = (nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++)
      vsigma[i] = 0.0;
  }

  if(v2rho2 != NULL){
    assert(info->provides & XC_PROVIDES_FXC);
    assert(v2rhosigma!=NULL && v2sigma2!=NULL);

    n = (nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++)
      v2rho2[i] = 0.0;

    n = (nspin == XC_UNPOLARIZED) ? 1 : 6;
    for(i=0; i<n; i++){
      v2rhosigma[i] = 0.0;
      v2sigma2[i]   = 0.0;
    }
  }

  /* check if density is larger than threshold */
  dens = rho[0];
  if(nspin == XC_POLARIZED) dens += rho[1];
  return (dens > MIN_DENS);
 
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
void XC(gga)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
	     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  assert(p!=NULL && p->info!=NULL);
  
  if(!XC(gga_input_init)(p->info, p->nspin, rho, zk, vrho, vsigma,
			 v2rho2, v2rhosigma, v2sigma2)) return;

  /* call functional */
  if(p->info->gga != NULL)
    p->info->gga(p, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

  if(p->mix != NULL)
    XC(mix_func)(p->mix, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

/* especializations */
inline void XC(gga_exc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
			FLOAT *zk)
{
  XC(gga)(p, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL);
}

inline void XC(gga_vxc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
			FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga)(p, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL);
}

inline void XC(gga_fxc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
			FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)(p, rho, sigma, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2);
}
