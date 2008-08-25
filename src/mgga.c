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
#include "funcs_mgga.c"

/* initialization */
int XC(mgga_init)(XC(mgga_type) *p, int functional, int nspin)
{
  int i;

  /* sanity check */
  assert(p != NULL);
  
  /* let us first find out if we know the functional */
  for(i=0; XC(mgga_known_funct)[i]!=NULL; i++){
    if(XC(mgga_known_funct)[i]->number == functional) break;
  }
  if(XC(mgga_known_funct)[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->params = NULL;
  p->info = XC(mgga_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}


void XC(mgga_end)(XC(mgga_type) *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);
}


void XC(mgga)(const XC(mgga_type) *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
	      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau,
	      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  FLOAT dens;
  int i, n;

  assert(p!=NULL && p->info!=NULL);
  
  /* sanity check */
  if(zk != NULL && !(p->info->provides & XC_PROVIDES_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc",
	    p->info->name);
    exit(1);
  }

  if(vrho != NULL && !(p->info->provides & XC_PROVIDES_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc",
	    p->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(p->info->provides & XC_PROVIDES_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc",
	    p->info->name);
    exit(1);
  }

  /* initialize output to zero */
  if(zk != NULL){
    assert(p->info->provides & XC_PROVIDES_EXC);
    *zk = 0.0;
  }

  if(vrho != NULL){
    assert(p->info->provides & XC_PROVIDES_VXC);
    assert(vsigma != NULL);

    for(i=0; i<p->nspin; i++){
      vrho[i] = 0.0;
      vtau[i] = 0.0;
    }

    n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++)
      vsigma[i] = 0.0;
  }

  if(v2rho2 != NULL){
    assert(p->info->provides & XC_PROVIDES_FXC);
    assert(v2rhosigma!=NULL && v2sigma2!=NULL && v2rhotau!=NULL && v2tausigma!=NULL && v2tau2!=NULL);

    n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
    for(i=0; i<n; i++){
      v2rho2[i] = 0.0;
      v2tau2[i] = 0.0;
    }

    n = (p->nspin == XC_UNPOLARIZED) ? 1 : 6;
    for(i=0; i<n; i++){
      v2rhosigma[i] = 0.0;
      v2tausigma[i] = 0.0;
      v2sigma2[i]   = 0.0;
    }
  }

  /* check if density is larger than threshold */
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  if(dens <= MIN_DENS) return;

  /* call functional */
  assert(p->info->mgga != NULL);
  p->info->mgga(p, rho, sigma, tau, zk, vrho, vsigma, vtau, 
		v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}

/* especializations */
inline void XC(mgga_exc)(const XC(mgga_type) *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau, 
			FLOAT *zk)
{
  XC(mgga)(p, rho, sigma, tau, zk, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void XC(mgga_vxc)(const XC(mgga_type) *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
			FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau)
{
  XC(mgga)(p, rho, sigma, tau, zk, vrho, vsigma, vtau, NULL, NULL, NULL, NULL, NULL, NULL);
}

inline void XC(mgga_fxc)(const XC(mgga_type) *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
			 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga)(p, rho, sigma, tau, NULL, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}


