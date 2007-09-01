#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

extern xc_func_info_type /* these are the LDA functionals that I know */
  func_info_lda_x,
  func_info_lda_c_wigner, 
  func_info_lda_c_rpa,
  func_info_lda_c_hl,
  func_info_lda_c_gl,
  func_info_lda_c_xalpha,
  func_info_lda_c_vwn,
  func_info_lda_c_vwn_rpa,
  func_info_lda_c_pz,
  func_info_lda_c_pz_mod,
  func_info_lda_c_ob_pz,
  func_info_lda_c_pw,
  func_info_lda_c_pw_mod,
  func_info_lda_c_ob_pw,
  func_info_lda_c_amgb;

const xc_func_info_type *lda_known_funct[] = {
  &func_info_lda_x,
  &func_info_lda_c_wigner,
  &func_info_lda_c_rpa,
  &func_info_lda_c_hl,
  &func_info_lda_c_gl,
  &func_info_lda_c_xalpha,
  &func_info_lda_c_vwn,
  &func_info_lda_c_vwn_rpa,
  &func_info_lda_c_pz,
  &func_info_lda_c_pz_mod,
  &func_info_lda_c_ob_pz,
  &func_info_lda_c_pw,
  &func_info_lda_c_pw_mod,
  &func_info_lda_c_ob_pw,
  &func_info_lda_c_amgb,
  NULL
};


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
void xc_lda(const xc_lda_type *p, const double *rho, double *exc, double *vxc, double *fxc, double *kxc)
{
  double dens;

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

/* get the lda functional */
void xc_lda_sp(const xc_lda_type *p, const float *rho, float *exc, float *vxc, float *fxc, float *kxc)
{
  double drho[2];
  double dexc;

  double * pexc = NULL;
  double * pvxc = NULL;
  double * pfxc = NULL;
  double * pkxc = NULL;

  int ii;
  const int nspin = p->nspin;

  drho[0] = rho[0];
  if(nspin > 1) drho[1] = rho[1];
  
  /* Allocate space for return values in double precision */
  if(exc) pexc = &dexc;
  if(vxc) pvxc = (double *) malloc(nspin * sizeof(double));
  if(fxc) pfxc = (double *) malloc(nspin * nspin * sizeof(double));
  if(kxc) pkxc = (double *) malloc(nspin * nspin * nspin * sizeof(double));

  /* Call the double precision version */
  xc_lda(p, drho, pexc, pvxc, pfxc, pkxc);

  /* Copy the result to the single precision return values */
  if(exc) exc[0] = dexc;

  if(vxc) {
    for(ii = 0; ii < nspin; ii++) vxc[ii] = pvxc[ii];
    free(pvxc);
  }

  if(fxc) {
    for(ii = 0; ii < nspin*nspin; ii++) fxc[ii] = pfxc[ii];
    free(pfxc);
  }

  if(kxc) {
    for(ii = 0; ii < nspin*nspin*nspin; ii++) kxc[ii] = pkxc[ii];
    free(pkxc);
  }

}


/* especializations */
void xc_lda_exc(const xc_lda_type *p, const double *rho, double *exc)
{
  xc_lda(p, rho, exc, NULL, NULL, NULL);
}

void xc_lda_vxc(const xc_lda_type *p, const double *rho, double *exc, double *vxc)
{
  xc_lda(p, rho, exc, vxc, NULL, NULL);
}

void xc_lda_fxc(const xc_lda_type *p, const double *rho, double *fxc)
{
  xc_lda(p, rho, NULL, NULL, fxc, NULL);
}

void xc_lda_kxc(const xc_lda_type *p, const double *rho, double *kxc)
{
  xc_lda(p, rho, NULL, NULL, NULL, kxc);
}

/* especialization in single precision */
void xc_lda_exc_sp(const xc_lda_type *p, const float *rho, float *exc)
{
  xc_lda_sp(p, rho, exc, NULL, NULL, NULL);
}

void xc_lda_vxc_sp(const xc_lda_type *p, const float *rho, float *exc, float *vxc)
{
  xc_lda_sp(p, rho, exc, vxc, NULL, NULL);
}

void xc_lda_fxc_sp(const xc_lda_type *p, const float *rho, float *fxc)
{
  xc_lda_sp(p, rho, NULL, NULL, fxc, NULL);
}

void xc_lda_kxc_sp(const xc_lda_type *p, const float *rho, float *kxc)
{
  xc_lda_sp(p, rho, NULL, NULL, NULL, kxc);
}


/* get the xc kernel through finite differences */
void xc_lda_fxc_fd(const xc_lda_type *p, const double *rho, double *fxc)
{
  static const double delta_rho = 1e-8;
  int i;

  for(i=0; i<p->nspin; i++){
    double rho2[2], e, vc1[2], vc2[2];
    int j, js;

    j  = (i+1) % 2;
    js = (i==0) ? 0 : 2;

    rho2[i] = rho[i] + delta_rho;
    rho2[j] = rho[j];
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

void xc_lda_kxc_fd(const xc_lda_type *p, const double *rho, double *kxc)
{
  /* Kxc, this is a third order tensor with respect to the densities */

  int i, j, k;
  const double delta_rho = 1e-4;

  for(i=0; i < p->nspin; i++){
    for(j=0; j < p->nspin; j++){
      for(k=0; k < p->nspin; k++){
    

	double rho2[2], e, vc1[2], vc2[2], vc3[2];
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
