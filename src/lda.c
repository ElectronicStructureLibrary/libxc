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
int XC(lda_init)(XC(lda_type) *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; XC(lda_known_funct)[i]!=NULL; i++){
    if(XC(lda_known_funct)[i]->number == functional) break;
  }
  assert(XC(lda_known_funct)[i] != NULL);
  if(XC(lda_known_funct)[i] == NULL) return -1; /* functional not found */
  
  /* initialize structure */
  p->params = NULL;
  p->info   = XC(lda_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  p->relativistic = 0;

  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}


/* termination */
void XC(lda_end)(XC(lda_type) *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);
}


/* get the lda functional */
void XC(lda)(const XC(lda_type) *p, const FLOAT *rho, 
	    FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
{
  FLOAT dens;

  assert(p!=NULL);
  
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

  { /* initialize output to zero */
    int i;

    if(zk != NULL) *zk = 0.0;

    if(vrho != NULL){
      for(i=0; i<p->nspin; i++)
    	vrho[i] = 0.0;
    }

    if(v2rho2 != NULL){
      int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
      for(i=0; i<n; i++)
	v2rho2[i] = 0.0;
    }

    if(v3rho3 != NULL){
      int n = (p->nspin == XC_UNPOLARIZED) ? 1 : 4;
      for(i=0; i<n; i++)
    	v3rho3[i] = 0.0;
    }
  }

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];

  if(dens <= MIN_DENS)
    return;

  assert(p->info!=NULL && p->info->lda!=NULL);

  /* call the LDA routines */
  p->info->lda(p, rho, zk, vrho, v2rho2, v3rho3);

  /* if necessary, call the finite difference routines */
  if(v2rho2 != NULL && !(p->info->provides & XC_PROVIDES_FXC))
    XC(lda_fxc_fd)(p, rho, v2rho2);
  
  if(v3rho3 != NULL && !(p->info->provides & XC_PROVIDES_KXC))
    XC(lda_kxc_fd)(p, rho, v3rho3);
}


/* especializations */
inline void XC(lda_exc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *zk)
{
  XC(lda)(p, rho, zk, NULL, NULL, NULL);
}

inline void XC(lda_vxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *zk, FLOAT *vrho)
{
  XC(lda)(p, rho, zk, vrho, NULL, NULL);
}

inline void XC(lda_fxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *v2rho2)
{
  XC(lda)(p, rho, NULL, NULL, v2rho2, NULL);
}

inline void XC(lda_kxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *v3rho3)
{
  XC(lda)(p, rho, NULL, NULL, NULL, v3rho3);
}


/* get the xc kernel through finite differences */
void XC(lda_fxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *fxc)
{
#ifdef SINGLE_PRECISION
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
    XC(lda_vxc)(p, rho2, &e, vc1);

    if(rho[i]<2.0*delta_rho){ /* we have to use a forward difference */
      XC(lda_vxc)(p, rho, &e, vc2);
	
      fxc[js] = (vc1[i] - vc2[i])/(delta_rho);
      if(p->nspin == XC_POLARIZED && i==0)
	fxc[1] = (vc1[j] - vc2[j])/(delta_rho);
	
    }else{                    /* centered difference (more precise)  */
      rho2[i] = rho[i] - delta_rho;
      XC(lda_vxc)(p, rho2, &e, vc2);
      
      fxc[js] = (vc1[i] - vc2[i])/(2.0*delta_rho);
      if(p->nspin == XC_POLARIZED && i==0)
	fxc[1] = (vc1[j] - vc2[j])/(2.0*delta_rho);
    }
    
  }
}


void XC(lda_kxc_fd)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *kxc)
{
  /* Kxc, this is a third order tensor with respect to the densities */

  int i, j, n;
  const FLOAT delta_rho = 1e-4;

  for(i=0; i<p->nspin; i++){
    FLOAT rho2[2], e, vc1[2], vc2[2], vc3[2];

    for(n=0; n<p->nspin; n++) rho2[n] = rho[n];
    XC(lda_vxc)(p, rho , &e, vc2);

    rho2[i] += delta_rho;
    XC(lda_vxc)(p, rho2, &e, vc1);
	
    rho2[i] -= 2.0*delta_rho;
    XC(lda_vxc)(p, rho2, &e, vc3);    
    
    for(j=0; j<p->nspin; j++)
      kxc[i*p->nspin + j] = (vc1[j] - 2.0*vc2[j] + vc3[j])/(delta_rho*delta_rho);
  }
}
