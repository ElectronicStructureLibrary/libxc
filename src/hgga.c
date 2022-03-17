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
xc_hgga_sanity_check(const xc_func_info_type *info, int order, xc_output_variables *out)
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

void xc_evaluate_hgga(const xc_func_type *func, int order, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau, const double *exx,
                xc_output_variables *out)
{
  xc_hgga_sanity_check(func->info, order, out);
  xc_output_variables_initialize(out, np, func->nspin);

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

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, lapl, tau, exx, out);
}
