/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

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
  p->mix_coef   = (double *) libxc_malloc(n_funcs*sizeof(double));
  p->func_aux   = (xc_func_type **) libxc_malloc(n_funcs*sizeof(xc_func_type *));

  for(ii=0; ii<n_funcs; ii++){
    p->mix_coef[ii] = mix_coef[ii];
    p->func_aux[ii] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
    xc_func_init (p->func_aux[ii], funcs_id[ii], p->nspin);
  }

  /* initialize variables */
  p->cam_omega = 0.0;
  p->cam_alpha = 0.0;
  p->cam_beta  = 0.0;
  p->nlc_b     = 0.0;
  p->nlc_C     = 0.0;
}

#ifdef HAVE_CUDA
__global__ static void add_to_mix_gpu(size_t np, double * dst, double coeff, double *src){
  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  if(ip < np) dst[ip] += coeff*src[ip];
}
#endif

static void add_to_mix(size_t np, double * dst, double coeff, double *src){
#ifndef HAVE_CUDA
  size_t ip;
  for(ip = 0; ip < np; ip++) dst[ip] += coeff*src[ip];
#else
  auto nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;
  add_to_mix_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>(np, dst, coeff, src);
#endif
}

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || (id) == XC_FAMILY_HYB_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA  || (id) == XC_FAMILY_HYB_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA  || (id) == XC_FAMILY_HYB_LDA ||  is_gga(id))
#define safe_free(pt) if(pt != NULL) libxc_free(pt)
#define sum_var(VAR) add_to_mix(np*dim->VAR, VAR, func->mix_coef[ii], VAR ## _);

void
xc_mix_func(const xc_func_type *func, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk, MGGA_OUT_PARAMS_NO_EXC(double *))
{
  const xc_func_type *aux;
  double *zk_;
  double *vrho_, *vsigma_, *vlapl_, *vtau_;
  double *v2rho2_, *v2rhosigma_, *v2rholapl_, *v2rhotau_, *v2sigma2_,
    *v2sigmalapl_, *v2sigmatau_, *v2lapl2_, *v2lapltau_,  *v2tau2_;
  double *v3rho3_, *v3rho2sigma_, *v3rho2lapl_, *v3rho2tau_, *v3rhosigma2_,
    *v3rhosigmalapl_, *v3rhosigmatau_, *v3rholapl2_, *v3rholapltau_,
    *v3rhotau2_, *v3sigma3_, *v3sigma2lapl_, *v3sigma2tau_,
    *v3sigmalapl2_, *v3sigmalapltau_, *v3sigmatau2_, *v3lapl3_,
    *v3lapl2tau_, *v3lapltau2_, *v3tau3_;
  double *v4rho4_, *v4rho3sigma_, *v4rho3lapl_, *v4rho3tau_, *v4rho2sigma2_,
    *v4rho2sigmalapl_, *v4rho2sigmatau_, *v4rho2lapl2_, *v4rho2lapltau_,
    *v4rho2tau2_, *v4rhosigma3_, *v4rhosigma2lapl_, *v4rhosigma2tau_,
    *v4rhosigmalapl2_, *v4rhosigmalapltau_, *v4rhosigmatau2_,
    *v4rholapl3_, *v4rholapl2tau_, *v4rholapltau2_, *v4rhotau3_,
    *v4sigma4_, *v4sigma3lapl_, *v4sigma3tau_, *v4sigma2lapl2_,
    *v4sigma2lapltau_, *v4sigma2tau2_, *v4sigmalapl3_, *v4sigmalapl2tau_,
    *v4sigmalapltau2_, *v4sigmatau3_, *v4lapl4_, *v4lapl3tau_,
    *v4lapl2tau2_, *v4lapltau3_, *v4tau4_;

  int ii;

  const xc_dimensions *dim = &(func->dim);

  /* Sanity check: have we claimed the highest possible derivatives?
     First, check for the lowest common derivative (also need to make
     sure the derivatives have been compiled in!)
  */
  int have_vxc = XC_FLAGS_I_HAVE_VXC;
  int have_fxc = XC_FLAGS_I_HAVE_FXC;
  int have_kxc = XC_FLAGS_I_HAVE_KXC;
  int have_lxc = XC_FLAGS_I_HAVE_LXC;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(! (aux->info->flags & XC_FLAGS_HAVE_VXC))
      have_vxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_FXC))
      have_fxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_KXC))
      have_kxc = 0;
    if(! (aux->info->flags & XC_FLAGS_HAVE_LXC))
      have_lxc = 0;
  }
  /* Then, for the actual checks */
  assert(have_lxc == (func->info->flags & XC_FLAGS_I_HAVE_LXC));
  assert(have_kxc == (func->info->flags & XC_FLAGS_I_HAVE_KXC));
  assert(have_fxc == (func->info->flags & XC_FLAGS_I_HAVE_FXC));
  assert(have_vxc == (func->info->flags & XC_FLAGS_I_HAVE_VXC));

  /* Sanity check: if component needs the Laplacian, then the mix
     must require it too */
  int need_laplacian = 0;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      need_laplacian = XC_FLAGS_NEEDS_LAPLACIAN;
  }
  assert((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) == need_laplacian);
  
  /* Check compatibility of the individual components */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    /* Sanity check: if component is GGA or meta-GGA, mix functional
       must also be GGA or meta-GGA */
    if(is_gga(aux->info->family))
      assert(is_gga(func->info->family));
    if(is_mgga(aux->info->family) && !is_mgga(func->info->family))
      assert(is_mgga(func->info->family));
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
  }

  /* prepare buffers that will hold the results from the individual functionals */
  zk_ = NULL;

  vrho_ = vsigma_ = vlapl_ = vtau_ = NULL;

  v2rho2_ = v2rhosigma_ = v2rholapl_ = v2rhotau_ = v2sigma2_ =
    v2sigmalapl_ = v2sigmatau_ = v2lapl2_ = v2lapltau_ =  v2tau2_ = NULL;

  v3rho3_ = v3rho2sigma_ = v3rho2lapl_ = v3rho2tau_ = v3rhosigma2_ =
    v3rhosigmalapl_ = v3rhosigmatau_ = v3rholapl2_ = v3rholapltau_ =
    v3rhotau2_ = v3sigma3_ = v3sigma2lapl_ = v3sigma2tau_ =
    v3sigmalapl2_ = v3sigmalapltau_ = v3sigmatau2_ = v3lapl3_ =
    v3lapl2tau_ = v3lapltau2_ = v3tau3_ = NULL;

  v4rho4_ = v4rho3sigma_ = v4rho3lapl_ = v4rho3tau_ = v4rho2sigma2_ =
    v4rho2sigmalapl_ = v4rho2sigmatau_ = v4rho2lapl2_ = v4rho2lapltau_ =
    v4rho2tau2_ = v4rhosigma3_ = v4rhosigma2lapl_ = v4rhosigma2tau_ =
    v4rhosigmalapl2_ = v4rhosigmalapltau_ = v4rhosigmatau2_ =
    v4rholapl3_ = v4rholapl2tau_ = v4rholapltau2_ = v4rhotau3_ =
    v4sigma4_ = v4sigma3lapl_ = v4sigma3tau_ = v4sigma2lapl2_ =
    v4sigma2lapltau_ = v4sigma2tau2_ = v4sigmalapl3_ = v4sigmalapl2tau_ =
    v4sigmalapltau2_ = v4sigmatau3_ = v4lapl4_ = v4lapl3tau_ =
    v4lapl2tau2_ = v4lapltau3_ = v4tau4_ = NULL;

  /* allocate buffers */
  if(zk != NULL)
    zk_ = (double *) libxc_malloc(sizeof(double)*np*dim->zk);

  if(vrho != NULL){
    vrho_ = (double *) libxc_malloc(sizeof(double)*np*dim->vrho);
    if(is_gga(func->info->family)){
      vsigma_ = (double *) libxc_malloc(sizeof(double)*np*dim->vsigma);
    }
    if(is_mgga(func->info->family)){
      if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
        vlapl_ = (double *) libxc_malloc(sizeof(double)*np*dim->vlapl);
      }
      vtau_  = (double *) libxc_malloc(sizeof(double)*np*dim->vtau);
    }
  }

  if(v2rho2 != NULL){
    v2rho2_ = (double *) libxc_malloc(sizeof(double)*np*dim->v2rho2);
    if(is_gga(func->info->family)){
      v2rhosigma_  = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhosigma);
      v2sigma2_    = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigma2);
    }
    if(is_mgga(func->info->family)){
      v2rholapl_   = (double *) libxc_malloc(sizeof(double)*np*dim->v2rholapl);
      v2rhotau_    = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhotau);
      v2sigmalapl_ = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmalapl);
      v2sigmatau_  = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmatau);
      v2lapl2_     = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapl2);
      v2lapltau_   = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapltau);
      v2tau2_      = (double *) libxc_malloc(sizeof(double)*np*dim->v2tau2);
    }
  }

  if(v3rho3 != NULL){
    v3rho3_      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho3);
    if(is_gga(func->info->family)){
      v3rho2sigma_ = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2sigma);
      v3rhosigma2_ = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigma2);
      v3sigma3_    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma3);
    }
    if(is_mgga(func->info->family)){
      v3rho2lapl_     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2lapl);
      v3rho2tau_      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2tau);
      v3rhosigmalapl_ = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmalapl);
      v3rhosigmatau_  = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmatau);
      v3rholapl2_     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapl2);
      v3rholapltau_   = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapltau);
      v3rhotau2_      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhotau2);
      v3sigma2lapl_   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2lapl);
      v3sigma2tau_    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2tau);
      v3sigmalapl2_   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapl2);
      v3sigmalapltau_ = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapltau);
      v3sigmatau2_    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmatau2);
      v3lapl3_        = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl3);
      v3lapl2tau_     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl2tau);
      v3lapltau2_     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapltau2);
      v3tau3_         = (double *) libxc_malloc(sizeof(double)*np*dim->v3tau3);
    }
  }
  if(v4rho4 != NULL){
    v4rho4_            = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho4);
    if(is_gga(func->info->family)){
      v4rho3sigma_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3sigma);
      v4rho2sigma2_      = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigma2);
      v4rhosigma3_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma3);
      v4sigma4_          = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma4);
    }
    if(is_mgga(func->info->family)){
      v4rho3lapl_        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3lapl);
      v4rho3tau_         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3tau);
      v4rho2sigmalapl_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmalapl);
      v4rho2sigmatau_    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmatau);
      v4rho2lapl2_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapl2);
      v4rho2lapltau_     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapltau);
      v4rho2tau2_        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2tau2);
      v4rhosigma2lapl_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2lapl);
      v4rhosigma2tau_    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2tau);
      v4rhosigmalapl2_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapl2);
      v4rhosigmalapltau_ = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapltau);
      v4rhosigmatau2_    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmatau2);
      v4rholapl3_        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl3);
      v4rholapl2tau_     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl2tau);
      v4rholapltau2_     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapltau2);
      v4rhotau3_         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhotau3);
      v4sigma3lapl_      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3lapl);
      v4sigma3tau_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3tau);
      v4sigma2lapl2_     = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapl2);
      v4sigma2lapltau_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapltau);
      v4sigma2tau2_      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2tau2);
      v4sigmalapl3_      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl3);
      v4sigmalapl2tau_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl2tau);
      v4sigmalapltau2_   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapltau2);
      v4sigmatau3_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmatau3);
      v4lapl4_           = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl4);
      v4lapl3tau_        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl3tau);
      v4lapl2tau2_       = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl2tau2);
      v4lapltau3_        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapltau3);
      v4tau4_            = (double *) libxc_malloc(sizeof(double)*np*dim->v4tau4);
    }
  }

  /* Proceed by computing the mix */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];

    /* Evaluate the functional */
    switch(aux->info->family){
    case XC_FAMILY_LDA:
      xc_lda(aux, np, rho, zk_, vrho_, v2rho2_, v3rho3_, v4rho4_);
      break;
    case XC_FAMILY_GGA:
      xc_gga(aux, np, rho, sigma, zk_, vrho_, vsigma_,
             v2rho2_, v2rhosigma_, v2sigma2_,
             v3rho3_, v3rho2sigma_, v3rhosigma2_, v3sigma3_,
             v4rho4_, v4rho3sigma_, v4rho2sigma2_, v4rhosigma3_, v4sigma4_);
      break;
    case XC_FAMILY_MGGA:
      xc_mgga(aux, np, rho, sigma, lapl, tau,
              zk_,
              vrho_, vsigma_, vlapl_, vtau_,
              v2rho2_, v2rhosigma_, v2rholapl_, v2rhotau_, v2sigma2_,
              v2sigmalapl_, v2sigmatau_, v2lapl2_, v2lapltau_,  v2tau2_,
              v3rho3_, v3rho2sigma_, v3rho2lapl_, v3rho2tau_, v3rhosigma2_,
              v3rhosigmalapl_, v3rhosigmatau_, v3rholapl2_, v3rholapltau_,
              v3rhotau2_, v3sigma3_, v3sigma2lapl_, v3sigma2tau_,
              v3sigmalapl2_, v3sigmalapltau_, v3sigmatau2_, v3lapl3_,
              v3lapl2tau_, v3lapltau2_, v3tau3_,
              v4rho4_, v4rho3sigma_, v4rho3lapl_, v4rho3tau_, v4rho2sigma2_,
              v4rho2sigmalapl_, v4rho2sigmatau_, v4rho2lapl2_, v4rho2lapltau_,
              v4rho2tau2_, v4rhosigma3_, v4rhosigma2lapl_, v4rhosigma2tau_,
              v4rhosigmalapl2_, v4rhosigmalapltau_, v4rhosigmatau2_,
              v4rholapl3_, v4rholapl2tau_, v4rholapltau2_, v4rhotau3_,
              v4sigma4_, v4sigma3lapl_, v4sigma3tau_, v4sigma2lapl2_,
              v4sigma2lapltau_, v4sigma2tau2_, v4sigmalapl3_, v4sigmalapl2tau_,
              v4sigmalapltau2_, v4sigmatau3_, v4lapl4_, v4lapl3tau_,
              v4lapl2tau2_, v4lapltau3_, v4tau4_);
      break;
    }

    /* Do the mixing */
    if(zk != NULL) {
      sum_var(zk);
    }

 #ifndef XC_DONT_COMPILE_VXC
    if(vrho != NULL) {
      sum_var(vrho);

      if(is_gga(aux->info->family)) {
        sum_var(vsigma);
      }
      
      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(vlapl);
        }
        sum_var(vtau);
      }
    }

#ifndef XC_DONT_COMPILE_FXC
    if(v2rho2 != NULL){
      sum_var(v2rho2);

      if(is_gga(aux->info->family)) {
        sum_var(v2rhosigma);
        sum_var(v2sigma2);
      }

      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v2rholapl);
          sum_var(v2sigmalapl);
          sum_var(v2lapl2);
          sum_var(v2lapltau);
        }          
        sum_var(v2rhotau);
        sum_var(v2sigmatau);
        sum_var(v2tau2);
      }
    }

#ifndef XC_DONT_COMPILE_KXC    
    if(v3rho3 != NULL){
      sum_var(v3rho3);

      if(is_gga(aux->info->family)) {
        sum_var(v3rho2sigma);
        sum_var(v3rhosigma2);
        sum_var(v3sigma3);
      }
      
      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v3rho2lapl);
          sum_var(v3rhosigmalapl);
          sum_var(v3rholapl2);
          sum_var(v3rholapltau);
          sum_var(v3sigma2lapl);
          sum_var(v3sigmalapl2);
          sum_var(v3sigmalapltau);
          sum_var(v3lapl3);
          sum_var(v3lapl2tau);
          sum_var(v3lapltau2);
        }
        sum_var(v3rho2tau);
        sum_var(v3rhosigmatau);
        sum_var(v3rhotau2);
        sum_var(v3sigma2tau);
        sum_var(v3sigmatau2);
        sum_var(v3tau3);
      }
    }

#ifndef XC_DONT_COMPILE_LXC
    if(v4rho4 != NULL){
      sum_var(v4rho4);

      if(is_gga(aux->info->family)) {
        sum_var(v4rho3sigma);
        sum_var(v4rho2sigma2);
        sum_var(v4rhosigma3);
        sum_var(v4sigma4);
      }
      if(is_mgga(aux->info->family)) {
        if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
          sum_var(v4rho3lapl);
          sum_var(v4rho2sigmalapl);
          sum_var(v4rho2lapl2);
          sum_var(v4rho2lapltau);
          sum_var(v4rhosigma2lapl);
          sum_var(v4rhosigmalapl2);
          sum_var(v4rhosigmalapltau);
          sum_var(v4rholapl3);
          sum_var(v4rholapl2tau);
          sum_var(v4rholapltau2);
          sum_var(v4sigma3lapl);
          sum_var(v4sigma2lapl2);
          sum_var(v4sigma2lapltau);
          sum_var(v4sigmalapl3);
          sum_var(v4sigmalapl2tau);
          sum_var(v4sigmalapltau2);
          sum_var(v4lapl4);
          sum_var(v4lapl3tau);
          sum_var(v4lapl2tau2);
          sum_var(v4lapltau3);
        }
        sum_var(v4rho3tau);
        sum_var(v4rho2sigmatau);
        sum_var(v4rho2tau2);
        sum_var(v4rhosigma2tau);
        sum_var(v4rhosigmatau2);
        sum_var(v4rhotau3);
        sum_var(v4sigma3tau);
        sum_var(v4sigma2tau2);
        sum_var(v4sigmatau3);
        sum_var(v4tau4);
      }
    }
#endif
#endif
#endif
#endif
  } /* end functional loop */

  /* deallocate internal buffers */
  safe_free(zk_);
#ifndef XC_DONT_COMPILE_VXC
  safe_free(vrho_); safe_free(vsigma_); safe_free(vlapl_); safe_free(vtau_);
#ifndef XC_DONT_COMPILE_FXC
  safe_free(v2rho2_); safe_free(v2rhosigma_); safe_free(v2rholapl_); safe_free(v2rhotau_);
  safe_free(v2sigma2_); safe_free(v2sigmalapl_); safe_free(v2sigmatau_);
  safe_free(v2lapl2_); safe_free(v2lapltau_); safe_free(v2tau2_);
#ifndef XC_DONT_COMPILE_KXC
  safe_free(v3rho3_); safe_free(v3rho2sigma_); safe_free(v3rho2lapl_); safe_free(v3rho2tau_);
  safe_free(v3rhosigma2_); safe_free(v3rhosigmalapl_); safe_free(v3rhosigmatau_);
  safe_free(v3rholapl2_); safe_free(v3rholapltau_); safe_free(v3rhotau2_);
  safe_free(v3sigma3_); safe_free(v3sigma2lapl_); safe_free(v3sigma2tau_);
  safe_free(v3sigmalapl2_); safe_free(v3sigmalapltau_); safe_free(v3sigmatau2_);
  safe_free(v3lapl3_); safe_free(v3lapl2tau_); safe_free(v3lapltau2_); safe_free(v3tau3_);
#ifndef XC_DONT_COMPILE_LXC
  safe_free(v4rho4_); safe_free(v4rho3sigma_); safe_free(v4rho3lapl_); safe_free(v4rho3tau_);
  safe_free(v4rho2sigma2_); safe_free(v4rho2sigmalapl_); safe_free(v4rho2sigmatau_);
  safe_free(v4rho2lapl2_); safe_free(v4rho2lapltau_); safe_free(v4rho2tau2_);
  safe_free(v4rhosigma3_); safe_free(v4rhosigma2lapl_); safe_free(v4rhosigma2tau_);
  safe_free(v4rhosigmalapl2_); safe_free(v4rhosigmalapltau_); safe_free(v4rhosigmatau2_);
  safe_free(v4rholapl3_); safe_free(v4rholapl2tau_); safe_free(v4rholapltau2_); safe_free(v4rhotau3_);
  safe_free(v4sigma4_); safe_free(v4sigma3lapl_); safe_free(v4sigma3tau_); safe_free(v4sigma2lapl2_);
  safe_free(v4sigma2lapltau_); safe_free(v4sigma2tau2_); safe_free(v4sigmalapl3_); safe_free(v4sigmalapl2tau_);
  safe_free(v4sigmalapltau2_); safe_free(v4sigmatau3_); safe_free(v4lapl4_); safe_free(v4lapl3tau_);
  safe_free(v4lapl2tau2_); safe_free(v4lapltau3_); safe_free(v4tau4_);
#endif
#endif
#endif
#endif
}

