/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA ||  is_gga(id))

/* initializes the mixing */
void
xc_mix_init(xc_func_type *p, int n_funcs, const int *funcs_id, const double *mix_coef)
{
  const xc_func_type *aux;
  int have_vxc, have_fxc, have_kxc, have_lxc, need_laplacian;
  int ii;

  assert(p != NULL);
  assert(p->func_aux == NULL && p->mix_coef == NULL);

  /* allocate structures needed for mixed functional */
  p->n_func_aux = n_funcs;
  p->mix_coef   = (double *) libxc_malloc(n_funcs*sizeof(double));
  p->func_aux   = (xc_func_type **) libxc_malloc(n_funcs*sizeof(xc_func_type *));

  for(ii=0; ii<n_funcs; ii++){
    p->mix_coef[ii] = mix_coef[ii];
    p->func_aux[ii] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
    xc_func_init (p->func_aux[ii], funcs_id[ii], p->nspin);
  }

  /* initialize variables */
  p->hyb_number_terms = 0;
  p->hyb_type  = NULL;
  p->hyb_coeff = NULL;
  p->hyb_omega = NULL;

  p->nlc_b     = 0.0;
  p->nlc_C     = 0.0;

  /* Sanity check: have we claimed the highest possible derivatives?
     First, check for the lowest common derivative (also need to make
     sure the derivatives have been compiled in!)
  */
  have_vxc = XC_FLAGS_I_HAVE_VXC;
  have_fxc = XC_FLAGS_I_HAVE_FXC;
  have_kxc = XC_FLAGS_I_HAVE_KXC;
  have_lxc = XC_FLAGS_I_HAVE_LXC;
  for(ii=0; ii<p->n_func_aux; ii++){
    aux = p->func_aux[ii];
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
  assert(have_lxc == (p->info->flags & XC_FLAGS_I_HAVE_LXC));
  assert(have_kxc == (p->info->flags & XC_FLAGS_I_HAVE_KXC));
  assert(have_fxc == (p->info->flags & XC_FLAGS_I_HAVE_FXC));
  assert(have_vxc == (p->info->flags & XC_FLAGS_I_HAVE_VXC));

  /* Sanity check: if component needs the Laplacian, then the mix
     must require it too */
  need_laplacian = 0;
  for(ii=0; ii<p->n_func_aux; ii++){
    aux = p->func_aux[ii];
    if(aux->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      need_laplacian = XC_FLAGS_NEEDS_LAPLACIAN;
  }
  assert((p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) == need_laplacian);

  /* Check compatibility of the individual components */
  for(ii=0; ii<p->n_func_aux; ii++){
    aux = p->func_aux[ii];
    /* Sanity check: if component is GGA or meta-GGA, mix functional
       must also be GGA or meta-GGA */
    if(is_gga(aux->info->family))
      assert(is_gga(p->info->family));
    if(is_mgga(aux->info->family) && !is_mgga(p->info->family))
      assert(is_mgga(p->info->family));
    /* Sanity checks: if mix functional has higher derivatives, these
       must also exist in the individual components */
    if(p->info->flags & XC_FLAGS_HAVE_VXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_VXC);
    if(p->info->flags & XC_FLAGS_HAVE_FXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_FXC);
    if(p->info->flags & XC_FLAGS_HAVE_KXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_KXC);
    if(p->info->flags & XC_FLAGS_HAVE_LXC)
      assert(aux->info->flags & XC_FLAGS_HAVE_LXC);
  }

}

#ifdef HAVE_CUDA
__global__ static void add_to_mix_gpu(size_t np, double * dst, double coeff, const double *src){
  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  if(ip < np) dst[ip] += coeff*src[ip];
}
#endif

static void add_to_mix(size_t np, double * dst, double coeff, const double *src){
#ifndef HAVE_CUDA
  size_t ip;
  for(ip = 0; ip < np; ip++) dst[ip] += coeff*src[ip];
#else
  auto nblocks = np/CUDA_BLOCK_SIZE;
  if(np != nblocks*CUDA_BLOCK_SIZE) nblocks++;
  add_to_mix_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>(np, dst, coeff, src);
#endif
}

#define init_var(VAR) for(ii=0; ii<dim->VAR; ii++) x ## VAR[ii] = 0.0
#define  sum_var(VAR) add_to_mix(dim->VAR, &(VAR[ip*dim->VAR]), func->mix_coef[ii], x ## VAR)

void
xc_mix_func(const xc_func_type *func, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  const xc_dimensions *dim = &(func->dim);
  const xc_func_type *aux;
  size_t ip;
  int ii;

  /* static allocations to simplify code and help GPU performance */
  double *xzk, mzk[1];

  double *xvrho, mvrho[2];
  double xvsigma[3];
  double xvlapl[2], xvtau[2];

  double *xv2rho2, mv2rho2[3];
  double xv2rhosigma[2*3], xv2sigma2[6];
  double xv2rholapl[2*2], xv2rhotau[2*2], xv2sigmalapl[3*2], xv2sigmatau[3*2];
  double xv2lapl2[3], xv2lapltau[2*2], xv2tau2[3];

  double *xv3rho3, mv3rho3[4];
  double xv3rho2sigma[3*3], xv3rhosigma2[2*6], xv3sigma3[10];
  double xv3rho2lapl[3*2], xv3rho2tau[3*2], xv3rhosigmalapl[2*3*2], xv3rhosigmatau[2*3*2];
  double xv3rholapl2[2*3], xv3rholapltau[2*2*2], xv3rhotau2[2*3], xv3sigma2lapl[6*2];
  double xv3sigma2tau[6*2], xv3sigmalapl2[3*3], xv3sigmalapltau[3*2*2], xv3sigmatau2[3*3];
  double xv3lapl3[4], xv3lapl2tau[3*2], xv3lapltau2[2*3], xv3tau3[4];

  double *xv4rho4, mv4rho4[5];
  double xv4rho3sigma[4*3], xv4rho2sigma2[3*6], xv4rhosigma3[2*10], xv4sigma4[15];
  double xv4rho3lapl[4*2], xv4rho3tau[4*2], xv4rho2sigmalapl[3*3*2], xv4rho2sigmatau[3*3*2];
  double xv4rho2lapl2[3*3], xv4rho2lapltau[3*2*2], xv4rho2tau2[3*3], xv4rhosigma2lapl[3*6*2];
  double xv4rhosigma2tau[3*6*2], xv4rhosigmalapl2[2*3*3], xv4rhosigmalapltau[2*3*2*2], xv4rhosigmatau2[2*6*3];
  double xv4rholapl3[2*4], xv4rholapl2tau[2*3*2], xv4rholapltau2[2*2*3], xv4rhotau3[2*4];
  double xv4sigma3lapl[10*2], xv4sigma3tau[10*3], xv4sigma2lapl2[6*3], xv4sigma2lapltau[6*2*2];
  double xv4sigma2tau2[6*3], xv4sigmalapl3[3*4], xv4sigmalapl2tau[3*3*2], xv4sigmalapltau2[3*2*3];
  double xv4sigmatau3[3*4], xv4lapl4[5], xv4lapl3tau[4*2], xv4lapl2tau2[3*3];
  double xv4lapltau3[2*4], xv4tau4[5];

  /* make sure we have the right pointers, so that the worker
     functions know what to calculate */
  xzk     = (    zk != NULL) ?     mzk : NULL;
  xvrho   = (  vrho != NULL) ?   mvrho : NULL;
  xv2rho2 = (v2rho2 != NULL) ? mv2rho2 : NULL;
  xv3rho3 = (v3rho3 != NULL) ? mv3rho3 : NULL;
  xv4rho4 = (v4rho4 != NULL) ? mv4rho4 : NULL;
  
  /* Proceed by computing the mix */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];

    for(ip=0; ip<np; ip++){
      /* Evaluate the functional */
      switch(aux->info->family){
      case XC_FAMILY_LDA:
        xc_lda(aux, 1, &(rho[ip*dim->rho]),
               xzk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
        break;
      case XC_FAMILY_GGA:
        xc_gga(aux, 1, &(rho[ip*dim->rho]), &(sigma[ip*dim->sigma]),
               xzk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
        break;
      case XC_FAMILY_MGGA:
        xc_mgga(aux, 1, &(rho[ip*dim->rho]), &(sigma[ip*dim->sigma]), &(lapl[ip*dim->lapl]), &(tau[ip*dim->tau]),
                xzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, x));
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
    } /* end point loop */
  } /* end functional loop */
}

