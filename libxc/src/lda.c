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
#include "funcs_lda.c"


/* initialization */
int xc_lda_init(xc_lda_type *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; lda_known_funct[i]!=NULL; i++){
    if(lda_known_funct[i]->number == functional) break;
  }
  assert(lda_known_funct[i] != NULL);
  if(lda_known_funct[i] == NULL) return -1; /* functional not found */
  
  /* initialize structure */
  p->info = lda_known_funct[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  p->relativistic = 0;

  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}

/* termination */
void xc_lda_end(xc_lda_type *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);
}

/* get the lda functional */
void xc_lda(const xc_lda_type *p, const FLOAT *rho, 
	    FLOAT *exc, FLOAT *vxc, FLOAT *fxc, FLOAT *kxc)
{
  FLOAT dens;

  assert(p!=NULL);
  
  { /* initialize output to zero */
    int i;

    if(exc != NULL) *exc = 0.0;

    if(vxc != NULL){
      for(i=0; i<p->nspin; i++)
    	vxc[i] = 0.0;
    }

    if(fxc != NULL){
      int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
      for(i=0; i<n; i++)
	fxc[i] = 0.0;
    }

    if(kxc != NULL){
      int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 4;
      for(i=0; i<n; i++)
    	kxc[i] = 0.0;
    }
  }

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];

  if(dens <= MIN_DENS)
    return;

  assert(p->info!=NULL && p->info->lda!=NULL);
  if((exc != NULL || vxc !=NULL) || 
     (fxc != NULL && (p->info->provides & XC_PROVIDES_FXC)))
    p->info->lda(p, rho, exc, vxc, fxc);

  if(fxc != NULL && !(p->info->provides & XC_PROVIDES_FXC))
    xc_lda_fxc_fd(p, rho, fxc);

  if(kxc != NULL && !(p->info->provides & XC_PROVIDES_KXC))
    xc_lda_kxc_fd(p, rho, kxc);
}


/* especializations */
inline void xc_lda_exc(const xc_lda_type *p, const FLOAT *rho, FLOAT *exc)
{
  xc_lda(p, rho, exc, NULL, NULL, NULL);
}

inline void xc_lda_vxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *exc, FLOAT *vxc)
{
  xc_lda(p, rho, exc, vxc, NULL, NULL);
}

inline void xc_lda_fxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *fxc)
{
  xc_lda(p, rho, NULL, NULL, fxc, NULL);
}

inline void xc_lda_kxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *kxc)
{
  xc_lda(p, rho, NULL, NULL, NULL, kxc);
}


/* get the xc kernel through finite differences */
void xc_lda_fxc_fd(const xc_lda_type *p, const FLOAT *rho, FLOAT *fxc)
{
#if _SINGLE_PREC_
  static const FLOAT delta_rho = 1e-8;
#else
  static const FLOAT delta_rho = 1e-6;
#endif
  int i;

  for(i=0; i<p->nspin; i++){
    FLOAT rho2[2], e, vc1[2], vc2[2];
    int j, js;

    j  = (i+1) % 2;
    js = (i==0) ? 0 : 2;

    rho2[i] = rho[i] + delta_rho;
    rho2[j] = (p->nspin == XC_POLARIZED) ? rho[j] : 0.0;
    xc_lda_vxc(p, rho2, &e, vc1);

    if(rho[i]<2.0*delta_rho){ /* we have to use a forward difference */
      xc_lda_vxc(p, rho, &e, vc2);
	
      fxc[js] = (vc1[i] - vc2[i])/(delta_rho);
      if(p->nspin == XC_POLARIZED && i==0)
	fxc[1] = (vc1[j] - vc2[j])/(delta_rho);
	
    }else{                    /* centered difference (more precise)  */
      rho2[i] = rho[i] - delta_rho;
      xc_lda_vxc(p, rho2, &e, vc2);
      
      fxc[js] = (vc1[i] - vc2[i])/(2.0*delta_rho);
      if(p->nspin == XC_POLARIZED && i==0)
	fxc[1] = (vc1[j] - vc2[j])/(2.0*delta_rho);
    }
    
  }
}


/* WANRNING - get rid of this by using new definition of output variable kxc */
#define ___(i, j, k) [2*(2*i + j) + k] 

void xc_lda_kxc_fd(const xc_lda_type *p, const FLOAT *rho, FLOAT *kxc)
{
  /* Kxc, this is a third order tensor with respect to the densities */

  int i, j, k;
  const FLOAT delta_rho = 1e-4;

  for(i=0; i < p->nspin; i++){
    for(j=0; j < p->nspin; j++){
      for(k=0; k < p->nspin; k++){
    

	FLOAT rho2[2], e, vc1[2], vc2[2], vc3[2];
	int n, der_dir, func_dir;

	/* This is a bit tricky, we have to calculate a third mixed
	   partial derivative, one is calculated analitically (Vxc)
	   and the other two numerically. */
	   
	/* Here we reorder the indexes so that the numerical
	   derivative is taken with respect to the same variable */

	if(i!=j) {
	  if(j==k){
	    func_dir = i;
	    der_dir  = j;
	  }else{
	    func_dir = j;
	    der_dir  = i;
	  }
	}else{
	  func_dir = k;
	  der_dir  = j;
	}

	for(n=0; n< p->nspin; n++) rho2[n] = rho[n];

	xc_lda_vxc(p, rho , &e, vc2);

	rho2[der_dir] += delta_rho;
	xc_lda_vxc(p, rho2, &e, vc1);
	
	rho2[der_dir] -= 2.0*delta_rho;
	xc_lda_vxc(p, rho2, &e, vc3);

	kxc ___(i, j, k) = (vc1[func_dir] - 2.0*vc2[func_dir] + vc3[func_dir])/(delta_rho*delta_rho);
	
      }
    }
  }

}
