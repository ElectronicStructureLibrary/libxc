/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola

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
  p->cam_omega = 0.0;
  p->cam_alpha = 0.0;
  p->cam_beta  = 0.0;
  p->nlc_b     = 0.0;
  p->nlc_C     = 0.0;
}

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA  || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA  || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))
#define safe_free(pt) if(pt != NULL) free(pt)

void
xc_mix_func(const xc_func_type *func, int np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk, MGGA_OUT_PARAMS_NO_EXC(double *))
{
  const xc_func_type *aux;
  double *zk_;
  double *vrho_, *vsigma_, *vlapl_, *vtau_;
  double *v2rho2_, *v2rhosigma_, *v2rholapl_, *v2rhotau_;
  double *v2sigma2_, *v2sigmalapl_, *v2sigmatau_;
  double *v2lapl2_, *v2lapltau_;
  double *v2tau2_;
  double *v3rho3_, *v3rho2sigma_, *v3rho2lapl_, *v3rho2tau_;
  double *v3rhosigma2_, *v3rhosigmalapl_, *v3rhosigmatau_;
  double *v3rholapl2_, *v3rholapltau_;
  double *v3rhotau2_;
  double *v3sigma3_, *v3sigma2lapl_, *v3sigma2tau_;
  double *v3sigmalapl2_, *v3sigmalapltau_;
  double *v3sigmatau2_;
  double *v3lapl3_, *v3lapl2tau_;
  double *v3lapltau2_;
  double *v3tau3_;

  int ip, ii;

  const xc_dimensions *dim = &(func->dim);

  /* Sanity check: have we claimed the highest possible derivatives?
     First, check for the lowest common derivative (also need to make
     sure the derivatives have been compiled in!)
  */
  int have_vxc = XC_FLAGS_I_HAVE_VXC;
  int have_fxc = XC_FLAGS_I_HAVE_FXC;
  int have_kxc = XC_FLAGS_I_HAVE_KXC;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(! (aux->info->flags & XC_FLAGS_HAVE_VXC))
      have_vxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_FXC))
      have_fxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_KXC))
      have_kxc = 0;
  }
  /* Then, for the actual checks */
  assert(have_kxc == (func->info->flags & XC_FLAGS_HAVE_KXC));
  assert(have_fxc == (func->info->flags & XC_FLAGS_HAVE_FXC));
  assert(have_vxc == (func->info->flags & XC_FLAGS_HAVE_VXC));

  /* prepare buffers that will hold the results from the individual functionals */
  zk_ = NULL;
  vrho_ = vsigma_ = vlapl_ = vtau_ = NULL;
  v2rho2_ = v2rhosigma_ = v2rholapl_ = v2rhotau_ = NULL;
  v2sigma2_ =  v2sigmalapl_ = v2sigmatau_ = NULL;
  v2lapl2_ = v2lapltau_ = NULL;
  v2tau2_ = NULL;
  v3rho3_ = v3rho2sigma_ = v3rho2lapl_ = v3rho2tau_ = NULL;
  v3rhosigma2_ = v3rhosigmalapl_ = v3rhosigmatau_ = NULL;
  v3rholapl2_ = v3rholapltau_ = NULL;
  v3rhotau2_ = NULL;
  v3sigma3_ = v3sigma2lapl_ = v3sigma2tau_ = NULL;
  v3sigmalapl2_ = v3sigmalapltau_ = NULL;
  v3sigmatau2_ = NULL;
  v3lapl3_ = v3lapl2tau_ = NULL;
  v3lapltau2_ = NULL;
  v3tau3_ = NULL;

  if(zk != NULL)
    zk_ = (double *) malloc(sizeof(double)*np*dim->zk);

  if(vrho != NULL){
    vrho_ = (double *) malloc(sizeof(double)*np*dim->vrho);
    if(is_gga(func->info->family)){
      vsigma_ = (double *) malloc(sizeof(double)*np*dim->vsigma);
    }
    if(is_mgga(func->info->family)){
      /* At the moment we always allocate the derivatives involving
         the laplacian, as some parts of Libxc do not take into
         account the XC_FLAGS_NEEDS_LAPLACIAN flag.
         if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){ */
      vlapl_ = (double *) malloc(sizeof(double)*np*dim->vlapl);
      /* } */
      vtau_  = (double *) malloc(sizeof(double)*np*dim->vtau);
    }
  }

  if(v2rho2 != NULL){
    v2rho2_ = (double *) malloc(sizeof(double)*np*dim->v2rho2);
    if(is_gga(func->info->family)){
      v2rhosigma_  = (double *) malloc(sizeof(double)*np*dim->v2rhosigma);
      v2sigma2_    = (double *) malloc(sizeof(double)*np*dim->v2sigma2);
    }
    if(is_mgga(func->info->family)){
      v2rholapl_   = (double *) malloc(sizeof(double)*np*dim->v2rholapl);
      v2rhotau_    = (double *) malloc(sizeof(double)*np*dim->v2rhotau);
      v2sigmalapl_ = (double *) malloc(sizeof(double)*np*dim->v2sigmalapl);
      v2sigmatau_  = (double *) malloc(sizeof(double)*np*dim->v2sigmatau);
      v2lapl2_     = (double *) malloc(sizeof(double)*np*dim->v2lapl2);
      v2lapltau_   = (double *) malloc(sizeof(double)*np*dim->v2lapltau);
      v2tau2_      = (double *) malloc(sizeof(double)*np*dim->v2tau2);
    }
  }

  if(v3rho3 != NULL){
    v3rho3_      = (double *) malloc(sizeof(double)*np*dim->v3rho3);
    if(is_gga(func->info->family)){
      v3rho2sigma_ = (double *) malloc(sizeof(double)*np*dim->v3rho2sigma);
      v3rhosigma2_ = (double *) malloc(sizeof(double)*np*dim->v3rhosigma2);
      v3sigma3_    = (double *) malloc(sizeof(double)*np*dim->v3sigma3);
    }
    if(is_mgga(func->info->family)){
      v3rho2lapl_     = (double *) malloc(sizeof(double)*np*dim->v3rho2lapl);
      v3rho2tau_      = (double *) malloc(sizeof(double)*np*dim->v3rho2tau);
      v3rhosigmalapl_ = (double *) malloc(sizeof(double)*np*dim->v3rhosigmalapl);
      v3rhosigmatau_  = (double *) malloc(sizeof(double)*np*dim->v3rhosigmatau);
      v3rholapl2_     = (double *) malloc(sizeof(double)*np*dim->v3rholapl2);
      v3rholapltau_   = (double *) malloc(sizeof(double)*np*dim->v3rholapltau);
      v3rhotau2_      = (double *) malloc(sizeof(double)*np*dim->v3rhotau2);
      v3sigma2lapl_   = (double *) malloc(sizeof(double)*np*dim->v3sigma2lapl);
      v3sigma2tau_    = (double *) malloc(sizeof(double)*np*dim->v3sigma2tau);
      v3sigmalapl2_   = (double *) malloc(sizeof(double)*np*dim->v3sigmalapl2);
      v3sigmalapltau_ = (double *) malloc(sizeof(double)*np*dim->v3sigmalapltau);;
      v3sigmatau2_    = (double *) malloc(sizeof(double)*np*dim->v3sigmatau2);
      v3lapl3_        = (double *) malloc(sizeof(double)*np*dim->v3lapl3);
      v3lapl2tau_     = (double *) malloc(sizeof(double)*np*dim->v3lapl2tau);
      v3lapltau2_     = (double *) malloc(sizeof(double)*np*dim->v3lapltau2);
      v3tau3_         = (double *) malloc(sizeof(double)*np*dim->v3tau3);
    }
  }

  /* Proceed by computing the mix */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    /* Sanity check: if component is GGA or meta-GGA, mix functional
       must also be GGA or meta-GGA */
    if(is_gga(aux->info->family))
      assert(is_gga(func->info->family));
    if(is_mgga(aux->info->family) && !is_mgga(func->info->family))
      assert(is_mgga(func->info->family));
    /* Sanity check: if component needs the Laplacian, then the mix
       must require it too */
    if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      assert(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN);
    /* Sanity checks: if mix functional has higher derivatives, these
       must also exist in the individual components */
    if(func->info->flags & XC_FLAGS_HAVE_VXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_VXC);
    if(func->info->flags & XC_FLAGS_HAVE_FXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_FXC);
    if(func->info->flags & XC_FLAGS_HAVE_KXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_KXC);
    if(func->info->flags & XC_FLAGS_HAVE_LXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_LXC);

    /* Evaluate the functional */
    switch(aux->info->family){
    case XC_FAMILY_LDA:
      xc_lda(aux, np, rho, zk_, vrho_, v2rho2_, v3rho3_);
      break;
    case XC_FAMILY_GGA:
      xc_gga(aux, np, rho, sigma, zk_, vrho_, vsigma_,
             v2rho2_, v2rhosigma_, v2sigma2_,
             v3rho3_, v3rho2sigma_, v3rhosigma2_, v3sigma3_);
      break;
    case XC_FAMILY_MGGA:
      xc_mgga(aux, np, rho, sigma, lapl, tau,
              zk_,
              vrho_, vsigma_, vlapl_, vtau_,
              v2rho2_, v2rhosigma_, v2rholapl_, v2rhotau_,
              v2sigma2_, v2sigmalapl_, v2sigmatau_,
              v2lapl2_, v2lapltau_,
              v2tau2_,
              v3rho3_, v3rho2sigma_, v3rho2lapl_, v3rho2tau_,
              v3rhosigma2_, v3rhosigmalapl_, v3rhosigmatau_,
              v3rholapl2_, v3rholapltau_,
              v3rhotau2_,
              v3sigma3_, v3sigma2lapl_, v3sigma2tau_,
              v3sigmalapl2_, v3sigmalapltau_,
              v3sigmatau2_,
              v3lapl3_, v3lapl2tau_,
              v3lapltau2_,
              v3tau3_
              );
      break;
    }

    /* Do the mixing */
    if(zk != NULL) {
      for(ip = 0; ip < np*dim->zk; ip++)
        zk[ip] += func->mix_coef[ii] * zk_[ip];
    }

    if(vrho != NULL) {
      for(ip = 0; ip < np*dim->vrho; ip++)
        vrho[ip] += func->mix_coef[ii] * vrho_[ip];

      if(is_gga(aux->info->family)) {
        for(ip = 0; ip < np*dim->vsigma; ip++)
          vsigma[ip] += func->mix_coef[ii] * vsigma_[ip];
      }

      if(is_mgga(aux->info->family)) {
        for(ip = 0; ip < np*dim->vtau; ip++)
          vtau[ip] += func->mix_coef[ii] * vtau_[ip];
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          for(ip = 0; ip < np*dim->vlapl; ip++)
            vlapl[ip] += func->mix_coef[ii] * vlapl_[ip];
        }
      }
    }

    if(v2rho2 != NULL){
      for(ip = 0; ip < np*dim->v2rho2; ip++)
        v2rho2[ip] += func->mix_coef[ii] * v2rho2_[ip];

      if(is_gga(aux->info->family)) {
        for(ip = 0; ip < np*dim->v2rhosigma; ip++)
          v2rhosigma[ip] += func->mix_coef[ii] * v2rhosigma_[ip];
        for(ip = 0; ip < np*dim->v2sigma2; ip++)
          v2sigma2[ip] += func->mix_coef[ii] * v2sigma2_[ip];
      }

      if(is_mgga(aux->info->family)) {
        for(ip = 0; ip < np*dim->v2rhotau; ip++)
          v2rhotau[ip]    += func->mix_coef[ii] * v2rhotau_[ip];
        for(ip = 0; ip < np*dim->v2sigmatau; ip++)
          v2sigmatau[ip]  += func->mix_coef[ii] * v2sigmatau_[ip];
        for(ip = 0; ip < np*dim->v2tau2; ip++)
          v2tau2[ip]      += func->mix_coef[ii] * v2tau2_[ip];
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          for(ip = 0; ip < np*dim->v2rholapl; ip++)
            v2rholapl[ip]   += func->mix_coef[ii] * v2rholapl_[ip];
          for(ip = 0; ip < np*dim->v2lapl2; ip++)
            v2lapl2[ip]     += func->mix_coef[ii] * v2lapl2_[ip];
          for(ip = 0; ip < np*dim->v2sigmalapl; ip++)
            v2sigmalapl[ip] += func->mix_coef[ii] * v2sigmalapl_[ip];
          for(ip = 0; ip < np*dim->v2lapltau; ip++)
            v2lapltau[ip]   += func->mix_coef[ii] * v2lapltau_[ip];
        }
      }
    }

    if(v3rho3 != NULL){
      for(ip = 0; ip < np*dim->v3rho3; ip++)
        v3rho3[ip] += func->mix_coef[ii] * v3rho3_[ip];

      if(is_gga(aux->info->family)) {
        for(ip = 0; ip < np*dim->v3rho2sigma; ip++)
          v3rho2sigma[ip] += func->mix_coef[ii] * v3rho2sigma_[ip];
        for(ip = 0; ip < np*dim->v3rhosigma2; ip++)
          v3rhosigma2[ip] += func->mix_coef[ii] * v3rhosigma2_[ip];
        for(ip = 0; ip < np*dim->v3sigma3; ip++)
          v3sigma3[ip] += func->mix_coef[ii] * v3sigma3_[ip];
      }
      if(is_mgga(aux->info->family)) {
        for(ip = 0; ip < np*dim->v3rho2tau; ip++)
          v3rho2tau[ip] += func->mix_coef[ii] * v3rho2tau_[ip];
        for(ip = 0; ip < np*dim->v3sigmatau2; ip++)
          v3sigmatau2[ip] += func->mix_coef[ii] * v3sigmatau2_[ip];
        for(ip = 0; ip < np*dim->v3rhotau2; ip++)
          v3rhotau2[ip] += func->mix_coef[ii] * v3rhotau2_[ip];
        for(ip = 0; ip < np*dim->v3sigma2tau; ip++)
          v3sigma2tau[ip] += func->mix_coef[ii] * v3sigma2tau_[ip];
        for(ip = 0; ip < np*dim->v3sigmatau2; ip++)
          v3sigmatau2[ip] += func->mix_coef[ii] * v3sigmatau2_[ip];
        for(ip = 0; ip < np*dim->v3tau3; ip++)
          v3tau3[ip] += func->mix_coef[ii] * v3tau3_[ip];

        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          for(ip = 0; ip < np*dim->v3rho2lapl; ip++)
            v3rho2lapl[ip] += func->mix_coef[ii] * v3rho2lapl_[ip];
          for(ip = 0; ip < np*dim->v3sigmalapl2; ip++)
            v3sigmalapl2[ip] += func->mix_coef[ii] * v3sigmalapl2_[ip];
          for(ip = 0; ip < np*dim->v3rholapl2; ip++)
            v3rholapl2[ip] += func->mix_coef[ii] * v3rholapl2_[ip];
          for(ip = 0; ip < np*dim->v3rholapltau; ip++)
            v3rholapltau[ip] += func->mix_coef[ii] * v3rholapltau_[ip];
          for(ip = 0; ip < np*dim->v3sigma2lapl; ip++)
            v3sigma2lapl[ip] += func->mix_coef[ii] * v3sigma2lapl_[ip];
          for(ip = 0; ip < np*dim->v3sigmalapl2; ip++)
            v3sigmalapl2[ip] += func->mix_coef[ii] * v3sigmalapl2_[ip];
          for(ip = 0; ip < np*dim->v3sigmalapltau; ip++)
            v3sigmalapltau[ip] += func->mix_coef[ii] * v3sigmalapltau_[ip];
          for(ip = 0; ip < np*dim->v3lapl3; ip++)
            v3lapl3[ip] += func->mix_coef[ii] * v3lapl3_[ip];
          for(ip = 0; ip < np*dim->v3lapl2tau; ip++)
            v3lapl2tau[ip] += func->mix_coef[ii] * v3lapl2tau_[ip];
          for(ip = 0; ip < np*dim->v3lapltau2; ip++)
            v3lapltau2[ip] += func->mix_coef[ii] * v3lapltau2_[ip];
        }
      }
    }
  }

  /* deallocate internal buffers */
  safe_free(zk_);
  safe_free(vrho_); safe_free(vsigma_); safe_free(vlapl_); safe_free(vtau_);
  safe_free(v2rho2_); safe_free(v2rhosigma_); safe_free(v2rholapl_); safe_free(v2rhotau_);
  safe_free(v2sigma2_); safe_free(v2sigmalapl_); safe_free(v2sigmatau_);
  safe_free(v2lapl2_); safe_free(v2lapltau_);
  safe_free(v2tau2_);
  safe_free(v3rho3_); safe_free(v3rho2sigma_); safe_free(v3rho2lapl_); safe_free(v3rho2tau_);
  safe_free(v3rhosigma2_); safe_free(v3rhosigmalapl_); safe_free(v3rhosigmatau_);
  safe_free(v3rholapl2_); safe_free(v3rholapltau_);
  safe_free(v3rhotau2_);
  safe_free(v3sigma3_); safe_free(v3sigma2lapl_); safe_free(v3sigma2tau_);
  safe_free(v3sigmalapl2_); safe_free(v3sigmalapltau_);
  safe_free(v3sigmatau2_);
  safe_free(v3lapl3_); safe_free(v3lapl2tau_);
  safe_free(v3lapltau2_);
  safe_free(v3tau3_);
}
