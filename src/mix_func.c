/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

/* initializes the mixing */
void 
xc_mix_init(xc_func_type *p, int n_funcs, const int *funcs_id, const double *mix_coef)
{
  int ii;

  assert(p != NULL);
  assert(p->func_aux == NULL && p->mix_coef == NULL);

  /* allocate structures needed for */
  p->n_func_aux = n_funcs;
  p->mix_coef   = (double *) malloc(n_funcs*sizeof(double));
  p->func_aux   = (xc_func_type **) malloc(n_funcs*sizeof(xc_func_type *));

  for(ii=0; ii<n_funcs; ii++){
    p->mix_coef[ii] = mix_coef[ii];
    p->func_aux[ii] = (xc_func_type *) malloc(sizeof(xc_func_type));
    xc_func_init (p->func_aux[ii], funcs_id[ii], p->nspin);
  }

  /* initialize variables */
  p->cam_omega=0.0;
  p->cam_alpha=0.0;
  p->cam_beta=0.0;
  p->nlc_b=0.0;
  p->nlc_C=0.0;
}

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA  || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA  || is_gga(is))
#define safe_free(pt) if(pt != NULL) free(pt)

void
xc_mix_func(const xc_func_type *func, int np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
            double *v2rho2, double *v2sigma2, double *v2lapl2, double *v2tau2,
            double *v2rhosigma, double *v2rholapl, double *v2rhotau, 
            double *v2sigmalapl, double *v2sigmatau, double *v2lapltau)
{
  const xc_func_type *aux;
  double *zk_, *vrho_, *vsigma_, *vlapl_, *vtau_;
  double *v2rho2_, *v2sigma2_, *v2lapl2_, *v2tau2_;
  double *v2rhosigma_, *v2rholapl_, *v2rhotau_;
  double *v2sigmalapl_, *v2sigmatau_, *v2lapltau_;

  int ip, ii;

  /* prepare buffers that will hold the results from the individual functionals */
  zk_ = NULL;
  vrho_ = vsigma_ = vlapl_ = vtau_ = NULL;
  v2rho2_ = v2sigma2_ = v2lapl2_ = v2tau2_ = NULL;
  v2rhosigma_ = v2rholapl_ = v2rhotau_ = NULL;
  v2sigmalapl_ = v2sigmatau_ = v2lapltau_ = NULL;

  if(zk != NULL)
    zk_ = (double *) malloc(sizeof(double)*np*func->n_zk);

  if(vrho != NULL){
    vrho_ = (double *) malloc(sizeof(double)*np*func->n_vrho);
    if(is_gga(func->info->family)){
      vsigma_ = (double *) malloc(sizeof(double)*np*func->n_vsigma);
    }
    if(is_mgga(func->info->family)){
      vtau_  = (double *) malloc(sizeof(double)*np*func->n_vtau);
      if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
        vlapl_ = (double *) malloc(sizeof(double)*np*func->n_vlapl);
    }
  }

  if(v2rho2 != NULL){
    v2rho2_ = (double *) malloc(sizeof(double)*np*func->n_v2rho2);
    if(is_gga(func->info->family)){
      v2sigma2_    = (double *) malloc(sizeof(double)*np*func->n_v2sigma2);
      v2rhosigma_  = (double *) malloc(sizeof(double)*np*func->n_v2rhosigma);
    }
    if(is_mgga(func->info->family)){
      v2lapl2_     = (double *) malloc(sizeof(double)*np*func->n_v2lapl2);
      v2tau2_      = (double *) malloc(sizeof(double)*np*func->n_v2tau2);
      v2rholapl_   = (double *) malloc(sizeof(double)*np*func->n_v2rholapl);
      v2rhotau_    = (double *) malloc(sizeof(double)*np*func->n_v2rhotau);
      v2sigmalapl_ = (double *) malloc(sizeof(double)*np*func->n_v2sigmalapl);
      v2sigmatau_  = (double *) malloc(sizeof(double)*np*func->n_v2sigmatau);
      v2lapltau_   = (double *) malloc(sizeof(double)*np*func->n_v2lapltau);
    }
  }

  /* we now add the different components */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    switch(aux->info->family){
    case XC_FAMILY_LDA:
      xc_lda(aux, np, rho, zk_, vrho_, v2rho2_, NULL);
      break;
    case XC_FAMILY_GGA:
      xc_gga(aux, np, rho, sigma, zk_, vrho_, vsigma_, v2rho2_, v2rhosigma_, v2sigma2_, NULL, NULL, NULL, NULL);
      break;
    case XC_FAMILY_MGGA:
      xc_mgga(aux, np, rho, sigma, lapl, tau, zk_, vrho_, vsigma_, vlapl_, vtau_,
	       v2rho2_, v2sigma2_, v2lapl2_, v2tau2_, v2rhosigma_, v2rholapl_, v2rhotau_, 
	       v2sigmalapl_, v2sigmatau_, v2lapltau_);
      break;
    }

    if(zk != NULL)
      for(ip = 0; ip < np*func->n_zk; ip++)
        zk[ip] += func->mix_coef[ii] * zk_[ip];

    if(vrho != NULL){
      for(ip = 0; ip < np*func->n_vrho; ip++)
        vrho[ip] += func->mix_coef[ii] * vrho_[ip];

      if(is_gga(func->info->family) && is_gga(aux->info->family))
        for(ip = 0; ip < np*func->n_vsigma; ip++)
          vsigma[ip] += func->mix_coef[ii] * vsigma_[ip];

      if(is_mgga(func->info->family) && is_mgga(aux->info->family)){
        for(ip = 0; ip < np*func->n_vtau; ip++)
          vtau[ip] += func->mix_coef[ii] * vtau_[ip];
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          /* Check that mix has been properly defined */
          assert(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN);
          for(ip = 0; ip < np*func->n_vlapl; ip++)
            vlapl[ip] += func->mix_coef[ii] * vlapl_[ip];
        }
      }
    }

    if(v2rho2 != NULL){
      for(ip = 0; ip < np*func->n_v2rho2; ip++)
        v2rho2[ip] += func->mix_coef[ii] * v2rho2_[ip];

      if(is_gga(func->info->family) && is_gga(aux->info->family)){
        for(ip = 0; ip < np*func->n_v2rhosigma; ip++)
          v2rhosigma[ip] += func->mix_coef[ii] * v2rhosigma_[ip];
        for(ip = 0; ip < np*func->n_v2sigma2; ip++)
          v2sigma2[ip] += func->mix_coef[ii] * v2sigma2_[ip];
      }

      if(is_mgga(func->info->family) && is_mgga(aux->info->family)){
        for(ip = 0; ip < np*func->n_v2lapl2; ip++)
          v2lapl2[ip]     += func->mix_coef[ii] * v2lapl2_[ip];
        for(ip = 0; ip < np*func->n_v2tau2; ip++)
          v2tau2[ip]      += func->mix_coef[ii] * v2tau2_[ip];
        for(ip = 0; ip < np*func->n_v2rholapl; ip++)
          v2rholapl[ip]   += func->mix_coef[ii] * v2rholapl_[ip];
        for(ip = 0; ip < np*func->n_v2rhotau; ip++)
          v2rhotau[ip]    += func->mix_coef[ii] * v2rhotau_[ip];
        for(ip = 0; ip < np*func->n_v2sigmalapl; ip++)
          v2sigmalapl[ip] += func->mix_coef[ii] * v2sigmalapl_[ip];
        for(ip = 0; ip < np*func->n_v2sigmatau; ip++)
          v2sigmatau[ip]  += func->mix_coef[ii] * v2sigmatau_[ip];
        for(ip = 0; ip < np*func->n_v2lapltau; ip++)
          v2lapltau[ip]   += func->mix_coef[ii] * v2lapltau_[ip];
      }
    }
  }

  /* deallocate internal buffers */
  safe_free(zk_);
  safe_free(vrho_); safe_free(vsigma_); safe_free(vlapl_); safe_free(vtau_);
  safe_free(v2rho2_); safe_free(v2sigma2_); safe_free(v2lapl2_); safe_free(v2tau2_);
  safe_free(v2rhosigma_); safe_free(v2rholapl_); safe_free(v2rhotau_);
  safe_free(v2sigmalapl_); safe_free(v2sigmatau_); safe_free(v2lapltau_);
}
