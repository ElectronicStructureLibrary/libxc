/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_hgga.c"

/* macro to check is a buffer exists */
#define check_out_var(VAR) if(out->VAR == NULL){fprintf(stderr, "error: output variable, out->" #VAR ", is a null pointer\n"); exit(1);}

void
xc_hgga_sanity_check(const xc_func_info_type *info, int order, xc_hgga_out_params *out)
{
  /* sanity check */
  if(order < 0 || order > 4){
    fprintf(stderr, "Order of derivatives '%d' not implemented\n",
            order);
    exit(1);
  }
  
  /* sanity check */
  if(out->zk != NULL && !(info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
	    info->name);
    exit(1);
  }

  if(out->vrho != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_VXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
              info->name);
      exit(1);
    }
    check_out_var(vsigma);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(vlapl);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(vtau);
    }
    check_out_var(vexx);
  }

  if(out->v2rho2 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_FXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v2rhosigma); 
    check_out_var(v2sigma2);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v2rholapl);
      check_out_var(v2sigmalapl);
      check_out_var(v2lapl2);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v2rhotau);
      check_out_var(v2sigmatau);
      check_out_var(v2tau2);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v2lapltau);
    }
  }

  if(out->v3rho3 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_KXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v3rho2sigma);
    check_out_var(v3rhosigma2);
    check_out_var(v3sigma3);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v3rho2lapl);
      check_out_var(v3rhosigmalapl);
      check_out_var(v3rholapl2);
      check_out_var(v3sigma2lapl);
      check_out_var(v3sigmalapl2);
      check_out_var(v3lapl3);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v3rho2tau);
      check_out_var(v3rhosigmatau);
      check_out_var(v3rhotau2);
      check_out_var(v3sigma2tau);
      check_out_var(v3sigmatau2);
      check_out_var(v3tau3);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v3rholapltau);
      check_out_var(v3sigmalapltau);
      check_out_var(v3lapl2tau);
      check_out_var(v3lapltau2);
    }
  }

  if(out->v4rho4 != NULL){
    if(!(info->flags & XC_FLAGS_HAVE_LXC)){
      fprintf(stderr, "Functional '%s' does not provide an implementation of lxc\n",
              info->name);
      exit(1);
    }
    check_out_var(v4rho3sigma);
    check_out_var(v4rho2sigma2);
    check_out_var(v4rhosigma3);
    check_out_var(v4sigma4);
    if(info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      check_out_var(v4rho3lapl);
      check_out_var(v4rho2sigmalapl);
      check_out_var(v4rho2lapl2);
      check_out_var(v4rhosigma2lapl);
      check_out_var(v4rhosigmalapl2);
      check_out_var(v4rholapl3);
      check_out_var(v4sigma3lapl);
      check_out_var(v4sigma2lapl2);
      check_out_var(v4sigmalapl3);
      check_out_var(v4lapl4);
    }
    if(info->flags & XC_FLAGS_NEEDS_TAU){
      check_out_var(v4rho3tau);
      check_out_var(v4rho2sigmatau);
      check_out_var(v4rho2tau2);
      check_out_var(v4rhosigma2tau);
      check_out_var(v4rhosigmatau2);
      check_out_var(v4rhotau3);
      check_out_var(v4sigma3tau);
      check_out_var(v4sigma2tau2);
      check_out_var(v4sigmatau3);
      check_out_var(v4tau4);
    }
    if((info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (info->flags & XC_FLAGS_NEEDS_TAU)){
      check_out_var(v4rho2lapltau);
      check_out_var(v4rhosigmalapltau);
      check_out_var(v4rholapl2tau);
      check_out_var(v4rholapltau2);
      check_out_var(v4sigma2lapltau);
      check_out_var(v4sigmalapl2tau);
      check_out_var(v4sigmalapltau2);
      check_out_var(v4lapl3tau);
      check_out_var(v4lapl2tau2);
      check_out_var(v4lapltau3); 
    }
  }
}

void
xc_hgga_initalize(const xc_func_type *func, size_t np, xc_hgga_out_params *out)
{
  const xc_dimensions *dim = func->dim;

  /* initialize output to zero */
  if(out->zk != NULL)
    libxc_memset(out->zk, 0, dim->zk*np*sizeof(double));

  if(out->vrho != NULL){
    libxc_memset(out->vrho,   0, dim->vrho  *np*sizeof(double));
    libxc_memset(out->vsigma, 0, dim->vsigma*np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) {
      libxc_memset(out->vlapl,  0, dim->vlapl *np*sizeof(double));
    }
    if(func->info->flags & XC_FLAGS_NEEDS_TAU) {
      libxc_memset(out->vtau,   0, dim->vtau  *np*sizeof(double));
    }
  }

  if(out->v2rho2 != NULL){
    libxc_memset(out->v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    libxc_memset(out->v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
    libxc_memset(out->v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      libxc_memset(out->v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      libxc_memset(out->v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
    }
    
    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
      libxc_memset(out->v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
      libxc_memset(out->v2tau2,     0, dim->v2tau2     *np*sizeof(double));
    }
    
    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
    }
  }

  if(out->v3rho3 != NULL){
    libxc_memset(out->v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    libxc_memset(out->v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
    libxc_memset(out->v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
    libxc_memset(out->v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      libxc_memset(out->v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));
      libxc_memset(out->v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      libxc_memset(out->v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      libxc_memset(out->v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      libxc_memset(out->v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
    }

    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
      libxc_memset(out->v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
      libxc_memset(out->v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
      libxc_memset(out->v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
      libxc_memset(out->v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));
      libxc_memset(out->v3tau3,        0, dim->v3tau3       *np*sizeof(double));
    }

    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      libxc_memset(out->v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      libxc_memset(out->v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
      libxc_memset(out->v3lapltau2,     0, dim->v3lapltau2    *np*sizeof(double));
    }
  }

  if(out->v4rho4 != NULL){
    libxc_memset(out->v4rho4,         0, dim->v4rho4        *np*sizeof(double));
    libxc_memset(out->v4rho3sigma,    0, dim->v4rho3sigma   *np*sizeof(double));
    libxc_memset(out->v4rho2sigma2,   0, dim->v4rho2sigma2  *np*sizeof(double));
    libxc_memset(out->v4rhosigma3,    0, dim->v4rhosigma3   *np*sizeof(double));
    libxc_memset(out->v4sigma4,       0, dim->v4sigma4      *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      libxc_memset(out->v4rho3lapl,        0, dim->v4rho3lapl       *np*sizeof(double));
      libxc_memset(out->v4rho2sigmalapl,   0, dim->v4rho2sigmalapl  *np*sizeof(double));
      libxc_memset(out->v4rho2lapl2,       0, dim->v4rho2lapl2      *np*sizeof(double));
      libxc_memset(out->v4rhosigma2lapl,   0, dim->v4rhosigma2lapl  *np*sizeof(double));
      libxc_memset(out->v4rhosigmalapl2,   0, dim->v4rhosigmalapl2  *np*sizeof(double));
      libxc_memset(out->v4rholapl3,        0, dim->v4rholapl3       *np*sizeof(double));
      libxc_memset(out->v4sigma3lapl,      0, dim->v4sigma3lapl     *np*sizeof(double));
      libxc_memset(out->v4sigma2lapl2,     0, dim->v4sigma2lapl2    *np*sizeof(double));
      libxc_memset(out->v4sigmalapl3,      0, dim->v4sigmalapl3     *np*sizeof(double));
      libxc_memset(out->v4lapl4,           0, dim->v4lapl4          *np*sizeof(double));
    }

    if(func->info->flags & XC_FLAGS_NEEDS_TAU){
      libxc_memset(out->v4rho3tau,      0, dim->v4rho3tau     *np*sizeof(double));
      libxc_memset(out->v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      libxc_memset(out->v4rho2tau2,     0, dim->v4rho2tau2    *np*sizeof(double));
      libxc_memset(out->v4rhosigma2tau, 0, dim->v4rhosigma2tau*np*sizeof(double));
      libxc_memset(out->v4rhosigmatau2, 0, dim->v4rhosigmatau2*np*sizeof(double));
      libxc_memset(out->v4rhotau3,      0, dim->v4rhotau3     *np*sizeof(double));
      libxc_memset(out->v4sigma3tau,    0, dim->v4sigma3tau   *np*sizeof(double));
      libxc_memset(out->v4sigma2tau2,   0, dim->v4sigma2tau2  *np*sizeof(double));
      libxc_memset(out->v4sigmatau3,    0, dim->v4sigmatau3   *np*sizeof(double));
      libxc_memset(out->v4tau4,         0, dim->v4tau4        *np*sizeof(double));
    }

    if((func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN) && (func->info->flags & XC_FLAGS_NEEDS_TAU)) {
      libxc_memset(out->v4rho2lapltau,     0, dim->v4rho2lapltau    *np*sizeof(double));
      libxc_memset(out->v4rhosigmalapltau, 0, dim->v4rhosigmalapltau*np*sizeof(double));
      libxc_memset(out->v4rholapl2tau,     0, dim->v4rholapl2tau    *np*sizeof(double));
      libxc_memset(out->v4rholapltau2,     0, dim->v4rholapltau2    *np*sizeof(double));
      libxc_memset(out->v4sigma2lapltau,   0, dim->v4sigma2lapltau  *np*sizeof(double));
      libxc_memset(out->v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      libxc_memset(out->v4sigmalapltau2,   0, dim->v4sigmalapltau2  *np*sizeof(double));
      libxc_memset(out->v4lapl3tau,        0, dim->v4lapl3tau       *np*sizeof(double));
      libxc_memset(out->v4lapl2tau2,       0, dim->v4lapl2tau2      *np*sizeof(double));
      libxc_memset(out->v4lapltau3,        0, dim->v4lapltau3       *np*sizeof(double));
    }
  }

}

void xc_hgga_new(const xc_func_type *func, int order, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau, const double *exx,
                xc_hgga_out_params *out)
{
  xc_hgga_sanity_check(func->info, order, out);
  xc_hgga_initalize(func, np, out);

  /* call the hGGA routines */
  if(func->info->hgga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->hgga->unpol[order] != NULL)
        func->info->hgga->unpol[order](func, np, rho, sigma, lapl, tau, exx, out);
    }else{
      if(func->info->hgga->pol[order] != NULL)
        func->info->hgga->pol[order](func, np, rho, sigma, lapl, tau, exx, out);
    }
  }

  /* WARNING no mixed functionals for now */
}
