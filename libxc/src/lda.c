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
void xc_lda(xc_lda_type *p, double *rho, double *ec, double *vc, double *fc)
{
  double dens;

  assert(p!=NULL);
  
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];

  if(dens <= MIN_DENS){
    int i;
      
    *ec = 0.0;
    for(i=0; i<p->nspin; i++) vc[i] = 0.0;
    return;
  }
    
  assert(p->info!=NULL && p->info->lda!=NULL);
  p->info->lda(p, rho, ec, vc, fc);
}

/* get the xc kernel */
void xc_lda_fxc(xc_lda_type *p, double *rho, double *fxc)
{
  if(p->info->provides & XC_PROVIDES_FXC){
    double ec, vc[2];
    xc_lda(p, rho, &ec, vc, fxc);

  }else{ /* get fxc through a numerical derivative */
    int i, j;
    double delta_rho = 1e-8;

    for(i=0; i<p->nspin; i++){
      double rho2[2], e, vc1[2], vc2[2];

      j = (i+1) % 2;

      rho2[i] = rho[i] + delta_rho;
      rho2[j] = rho[j];
      xc_lda(p, rho2, &e, vc1, NULL);

      if(rho[i]<2.0*delta_rho){ /* we have to use a forward difference */
	xc_lda(p, rho, &e, vc2, NULL);
	
	fxc __(i, i) = (vc1[i] - vc2[i])/(delta_rho);
	if(p->nspin == XC_POLARIZED)
	  fxc __(i, j) = (vc1[j] - vc2[j])/(delta_rho);
	
      }else{                    /* centered difference (more precise)  */
	rho2[i] = rho[i] - delta_rho;
	xc_lda(p, rho2, &e, vc2, NULL);
	
	fxc __(i, i) = (vc1[i] - vc2[i])/(2.0*delta_rho);
	if(p->nspin == XC_POLARIZED)
	  fxc __(i, j) = (vc1[j] - vc2[j])/(2.0*delta_rho);
      }

    }
  }
}


void xc_lda_kxc(xc_lda_type *p, double *rho, double *kxc)
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

	xc_lda(p, rho , &e, vc2, NULL);

	rho2[der_dir] += delta_rho;
	xc_lda(p, rho2, &e, vc1, NULL);
	
	rho2[der_dir] -= 2.0*delta_rho;
	xc_lda(p, rho2, &e, vc3, NULL);

	kxc ___(i, j, k) = (vc1[func_dir] - 2.0*vc2[func_dir] + vc3[func_dir])/(delta_rho*delta_rho);
	
      }
    }
  }

}
