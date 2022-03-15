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
#define safe_free(pt) if(pt != NULL) libxc_free(pt)

/* Due to the fact that some functionals do not depend on tau or lapl,
   some variables may remain uninitialized when evaluating the ked or
   mgga functionals. Therefore, *all* variables should be explicitly
   initialized to zero here!
*/
void
xc_mgga_vars_allocate_all(int family, size_t np, const xc_dimensions *dim,
                     int do_zk, int do_vrho, int do_v2rho2, int do_v3rho3, int do_v4rho4,
                     double **zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double **, ))
{
  /* allocate buffers */
  if(do_zk){
    *zk = (double *) libxc_malloc(sizeof(double)*np*dim->zk);
    libxc_memset(*zk, 0, dim->zk*np*sizeof(double));
  }

#ifndef XC_DONT_COMPILE_VXC
  if(do_vrho){
    *vrho = (double *) libxc_malloc(sizeof(double)*np*dim->vrho);
    libxc_memset(*vrho,   0, dim->vrho  *np*sizeof(double));
    if(is_gga(family)){
      *vsigma = (double *) libxc_malloc(sizeof(double)*np*dim->vsigma);
      libxc_memset(*vsigma, 0, dim->vsigma*np*sizeof(double));
    }
    if(is_mgga(family)){
      *vlapl = (double *) libxc_malloc(sizeof(double)*np*dim->vlapl);
      libxc_memset(*vlapl,  0, dim->vlapl *np*sizeof(double));
      *vtau  = (double *) libxc_malloc(sizeof(double)*np*dim->vtau);
      libxc_memset(*vtau,   0, dim->vtau  *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_FXC
  if(do_v2rho2){
    *v2rho2 = (double *) libxc_malloc(sizeof(double)*np*dim->v2rho2);
     libxc_memset(*v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    if(is_gga(family)){
      *v2rhosigma  = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhosigma);
      libxc_memset(*v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
      *v2sigma2    = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigma2);
      libxc_memset(*v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v2rholapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v2rholapl);
      libxc_memset(*v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      *v2rhotau    = (double *) libxc_malloc(sizeof(double)*np*dim->v2rhotau);
      libxc_memset(*v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
      *v2sigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmalapl);
      libxc_memset(*v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      *v2sigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v2sigmatau);
      libxc_memset(*v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
      *v2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapl2);
      libxc_memset(*v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
      *v2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v2lapltau);
      libxc_memset(*v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
      *v2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v2tau2);
      libxc_memset(*v2tau2,     0, dim->v2tau2     *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_KXC
  if(do_v3rho3){
    *v3rho3      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho3);
    libxc_memset(*v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    if(is_gga(family)){
      *v3rho2sigma = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2sigma);
      libxc_memset(*v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
      *v3rhosigma2 = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigma2);
      libxc_memset(*v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
      *v3sigma3    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma3);
      libxc_memset(*v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v3rho2lapl     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2lapl);
      libxc_memset(*v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      *v3rho2tau      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rho2tau);
      libxc_memset(*v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
      *v3rhosigmalapl = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmalapl);
      libxc_memset(*v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));
      *v3rhosigmatau  = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhosigmatau);
      libxc_memset(*v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
      *v3rholapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapl2);
      libxc_memset(*v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      *v3rholapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v3rholapltau);
      libxc_memset(*v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      *v3rhotau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v3rhotau2);
      libxc_memset(*v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
      *v3sigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2lapl);
      libxc_memset(*v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      *v3sigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigma2tau);
      libxc_memset(*v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
      *v3sigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapl2);
      libxc_memset(*v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      *v3sigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmalapltau);
      libxc_memset(*v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      *v3sigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v3sigmatau2);
      libxc_memset(*v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));
      *v3lapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl3);
      libxc_memset(*v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
      *v3lapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapl2tau);
      libxc_memset(*v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
      *v3lapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v3lapltau2);
      libxc_memset(*v3lapltau2,     0, dim->v3lapltau2    *np*sizeof(double));
      *v3tau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v3tau3);
      libxc_memset(*v3tau3,        0, dim->v3tau3       *np*sizeof(double));
    }
  }

#ifndef XC_DONT_COMPILE_LXC
  if(do_v4rho4){
    *v4rho4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho4);
    libxc_memset(*v4rho4,         0, dim->v4rho4        *np*sizeof(double));
    if(is_gga(family)){
      *v4rho3sigma       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3sigma);
      libxc_memset(*v4rho3sigma,    0, dim->v4rho3sigma   *np*sizeof(double));
      *v4rho2sigma2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigma2);
      libxc_memset(*v4rho2sigma2,   0, dim->v4rho2sigma2  *np*sizeof(double));
      *v4rhosigma3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma3);
      libxc_memset(*v4rhosigma3,    0, dim->v4rhosigma3   *np*sizeof(double));
      *v4sigma4          = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma4);
      libxc_memset(*v4sigma4,       0, dim->v4sigma4      *np*sizeof(double));
    }
    if(is_mgga(family)){
      *v4rho3lapl        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3lapl);
      libxc_memset(*v4rho3lapl,        0, dim->v4rho3lapl       *np*sizeof(double));
      *v4rho3tau         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho3tau);
      libxc_memset(*v4rho3tau,      0, dim->v4rho3tau     *np*sizeof(double));
      *v4rho2sigmalapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmalapl);
      libxc_memset(*v4rho2sigmalapl,   0, dim->v4rho2sigmalapl  *np*sizeof(double));
      *v4rho2sigmatau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2sigmatau);
      libxc_memset(*v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      *v4rho2lapl2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapl2);
      libxc_memset(*v4rho2lapl2,       0, dim->v4rho2lapl2      *np*sizeof(double));
      *v4rho2lapltau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2lapltau);
      libxc_memset(*v4rho2lapltau,     0, dim->v4rho2lapltau    *np*sizeof(double));
      *v4rho2tau2        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rho2tau2);
      libxc_memset(*v4rho2tau2,     0, dim->v4rho2tau2    *np*sizeof(double));
      *v4rhosigma2lapl   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2lapl);
      libxc_memset(*v4rhosigma2lapl,   0, dim->v4rhosigma2lapl  *np*sizeof(double));
      *v4rhosigma2tau    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigma2tau);
      libxc_memset(*v4rho2sigmatau, 0, dim->v4rho2sigmatau*np*sizeof(double));
      *v4rhosigmalapl2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapl2);
      libxc_memset(*v4rhosigmalapl2,   0, dim->v4rhosigmalapl2  *np*sizeof(double));
      *v4rhosigmalapltau = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmalapltau);
      libxc_memset(*v4rhosigmalapltau, 0, dim->v4rhosigmalapltau*np*sizeof(double));
      *v4rhosigmatau2    = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhosigmatau2);
      libxc_memset(*v4rhosigmatau2, 0, dim->v4rhosigmatau2*np*sizeof(double));
      *v4rholapl3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl3);
      libxc_memset(*v4rholapl3,        0, dim->v4rholapl3       *np*sizeof(double));
      *v4rholapl2tau     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapl2tau);
      libxc_memset(*v4rholapl2tau,     0, dim->v4rholapl2tau    *np*sizeof(double));
      *v4rholapltau2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4rholapltau2);
      libxc_memset(*v4rholapltau2,     0, dim->v4rholapltau2    *np*sizeof(double));
      *v4rhotau3         = (double *) libxc_malloc(sizeof(double)*np*dim->v4rhotau3);
      libxc_memset(*v4rhotau3,      0, dim->v4rhotau3     *np*sizeof(double));
      *v4sigma3lapl      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3lapl);
      libxc_memset(*v4sigma3lapl,      0, dim->v4sigma3lapl     *np*sizeof(double));
      *v4sigma3tau       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma3tau);
      libxc_memset(*v4sigma3tau,    0, dim->v4sigma3tau   *np*sizeof(double));
      *v4sigma2lapl2     = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapl2);
      libxc_memset(*v4sigma2lapl2,     0, dim->v4sigma2lapl2    *np*sizeof(double));
      *v4sigma2lapltau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2lapltau);
      libxc_memset(*v4sigma2lapltau,   0, dim->v4sigma2lapltau  *np*sizeof(double));
      *v4sigma2tau2      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigma2tau2);
      libxc_memset(*v4sigma2tau2,   0, dim->v4sigma2tau2  *np*sizeof(double));
      *v4sigmalapl3      = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl3);
      libxc_memset(*v4sigmalapl3,      0, dim->v4sigmalapl3     *np*sizeof(double));
      *v4sigmalapl2tau   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapl2tau);
      libxc_memset(*v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      *v4sigmalapltau2   = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmalapltau2);
      libxc_memset(*v4sigmalapl2tau,   0, dim->v4sigmalapl2tau  *np*sizeof(double));
      *v4sigmatau3       = (double *) libxc_malloc(sizeof(double)*np*dim->v4sigmatau3);
      libxc_memset(*v4sigmatau3,    0, dim->v4sigmatau3   *np*sizeof(double));
      *v4lapl4           = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl4);
      libxc_memset(*v4lapl4,           0, dim->v4lapl4          *np*sizeof(double));
      *v4lapl3tau        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl3tau);
      libxc_memset(*v4lapl3tau,        0, dim->v4lapl3tau       *np*sizeof(double));
      *v4lapl2tau2       = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapl2tau2);
      libxc_memset(*v4lapl2tau2,       0, dim->v4lapl2tau2      *np*sizeof(double));
      *v4lapltau3        = (double *) libxc_malloc(sizeof(double)*np*dim->v4lapltau3);
      libxc_memset(*v4lapltau3,        0, dim->v4lapltau3       *np*sizeof(double));
      *v4tau4            = (double *) libxc_malloc(sizeof(double)*np*dim->v4tau4);
      libxc_memset(*v4tau4,         0, dim->v4tau4        *np*sizeof(double));
    }
  }
#endif
#endif
#endif
#endif
}

void
xc_mgga_vars_free_all(double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  /* deallocate internal buffers */
  safe_free(zk);
#ifndef XC_DONT_COMPILE_VXC
  safe_free(vrho);
  safe_free(vsigma);
  safe_free(vlapl);
  safe_free(vtau);

#ifndef XC_DONT_COMPILE_FXC
  safe_free(v2rho2); safe_free(v2rhosigma); safe_free(v2rholapl); safe_free(v2rhotau);
  safe_free(v2sigma2); safe_free(v2sigmalapl); safe_free(v2sigmatau);
  safe_free(v2lapl2); safe_free(v2lapltau); safe_free(v2tau2);

#ifndef XC_DONT_COMPILE_KXC
  safe_free(v3rho3); safe_free(v3rho2sigma); safe_free(v3rho2lapl); safe_free(v3rho2tau);
  safe_free(v3rhosigma2); safe_free(v3rhosigmalapl); safe_free(v3rhosigmatau);
  safe_free(v3rholapl2); safe_free(v3rholapltau); safe_free(v3rhotau2);
  safe_free(v3sigma3); safe_free(v3sigma2lapl); safe_free(v3sigma2tau);
  safe_free(v3sigmalapl2); safe_free(v3sigmalapltau); safe_free(v3sigmatau2);
  safe_free(v3lapl3); safe_free(v3lapl2tau); safe_free(v3lapltau2); safe_free(v3tau3);

#ifndef XC_DONT_COMPILE_LXC
  safe_free(v4rho4); safe_free(v4rho3sigma); safe_free(v4rho3lapl); safe_free(v4rho3tau);
  safe_free(v4rho2sigma2); safe_free(v4rho2sigmalapl); safe_free(v4rho2sigmatau);
  safe_free(v4rho2lapl2); safe_free(v4rho2lapltau); safe_free(v4rho2tau2);
  safe_free(v4rhosigma3); safe_free(v4rhosigma2lapl); safe_free(v4rhosigma2tau);
  safe_free(v4rhosigmalapl2); safe_free(v4rhosigmalapltau); safe_free(v4rhosigmatau2);
  safe_free(v4rholapl3); safe_free(v4rholapl2tau); safe_free(v4rholapltau2); safe_free(v4rhotau3);
  safe_free(v4sigma4); safe_free(v4sigma3lapl); safe_free(v4sigma3tau); safe_free(v4sigma2lapl2);
  safe_free(v4sigma2lapltau); safe_free(v4sigma2tau2); safe_free(v4sigmalapl3); safe_free(v4sigmalapl2tau);
  safe_free(v4sigmalapltau2); safe_free(v4sigmatau3); safe_free(v4lapl4); safe_free(v4lapl3tau);
  safe_free(v4lapl2tau2); safe_free(v4lapltau3); safe_free(v4tau4);
#endif
#endif
#endif
#endif
}

void
xc_mgga_evaluate_functional(const xc_func_type *func, size_t np,
                            const double *rho, const double *sigma, const double *lapl, const double *tau,
                            double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  double *mzk = NULL;

  if(func->info->flags & XC_FLAGS_HAVE_EXC)
    mzk = zk;
  
  /* Evaluate the functional */
  switch(func->info->family){
  case XC_FAMILY_LDA:
    xc_lda(func, np, rho,
           mzk LDA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_GGA:
    xc_gga(func, np, rho, sigma,
           mzk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  case XC_FAMILY_MGGA:
    xc_mgga(func, np, rho, sigma, lapl, tau,
            mzk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, ));
    break;
  }
}

void xc_mgga_evaluate_functional_new
(const xc_func_type *func, int order, size_t np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 xc_output_variables *out)
{
  /* Evaluate the functional */
  switch(func->info->family){
  case XC_FAMILY_LDA:
    xc_lda_new(func, order, np, rho, out);
    break;
  case XC_FAMILY_GGA:
    xc_gga_new(func, order, np, rho, sigma, out);
    break;
  case XC_FAMILY_MGGA:
    xc_mgga_new(func, order, np, rho, sigma, lapl, tau, out);
    break;
  }
}

/* initializes the mixing */
void
xc_deorbitalize_init(xc_func_type *p, int mgga_id, int ked_id)
{
  assert(p != NULL && p->func_aux == NULL);

  /* allocate structures needed for */
  p->n_func_aux = 2;
  p->func_aux   = (xc_func_type **) libxc_malloc(2*sizeof(xc_func_type *));

  p->func_aux[0] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));
  p->func_aux[1] = (xc_func_type *) libxc_malloc(sizeof(xc_func_type));

  xc_func_init (p->func_aux[0], mgga_id, p->nspin);
  xc_func_init (p->func_aux[1], ked_id,  p->nspin);
}

#define VAR(var, ip, index)        var[ip*p->dim->var + index]
static void
deorbitalized_functional_work
(const xc_func_type *p, size_t np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 xc_output_variables *out)
{
  double *mrho, *msigma, *mlapl, *mtau;
  xc_output_variables *mgga, *ked1, *ked2;
  size_t ip;
  
  int ii, max_order;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  max_order = -1;
  for(ii=0; ii <= XC_MAXIMUM_ORDER; ii++){
    if(orders[ii]) max_order = ii;
  }
  /* we actually need all orders <= max_order, so we change orders */
  for(ii=0; ii <= max_order; ii++)
    orders[ii] = 1;
  
  /* 
     prepare buffers that will hold the results from the individual functionals 
     we have to declare them as XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU
     as the code accesses these values (even if they are zero)
  */
  mgga = xc_allocate_output_variables
    (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);
                                         
  ked1 = xc_allocate_output_variables
    (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);

  mtau   = (double *) libxc_malloc(np*p->dim->tau*sizeof(double));
  if(p->nspin == XC_POLARIZED){
    mrho   = (double *) libxc_malloc(np*p->dim->rho*sizeof(double));
    msigma = (double *) libxc_malloc(np*p->dim->sigma*sizeof(double));
    mlapl  = (double *) libxc_malloc(np*p->dim->lapl*sizeof(double));

    ked2 = xc_allocate_output_variables
      (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);
  }
  
  /* evaluate the kinetic energy functional */
  
  if(p->nspin == XC_UNPOLARIZED){
    xc_mgga_evaluate_functional_new
      (p->func_aux[1], max_order, np, rho, sigma, lapl, NULL, ked1);
  }else{
    for(ip=0; ip<np; ip++){
      mrho  [2*ip] = rho  [2*ip]; mrho  [2*ip+1] = 0.0;
      msigma[3*ip] = sigma[3*ip]; msigma[3*ip+1] = 0.0; msigma[3*ip+2] = 0.0;
      mlapl [2*ip] = lapl [2*ip]; mlapl [2*ip+1] = 0.0;
    }
    xc_mgga_evaluate_functional_new
      (p->func_aux[1], max_order, np, mrho, msigma, mlapl, NULL, ked1);

    for(ip=0; ip<np; ip++){
      mrho  [2*ip] = rho  [2*ip + 1];
      msigma[3*ip] = sigma[3*ip + 2];
      mlapl [2*ip] = lapl [2*ip + 1];
    }
    xc_mgga_evaluate_functional_new
      (p->func_aux[1], max_order, np, mrho, msigma, mlapl, NULL, ked2);
  }
  
  /* now evaluate the mgga functional */
  if(p->nspin == XC_UNPOLARIZED){
    for(ip=0; ip<np; ip++){
      mtau[ip] = rho[ip]*ked1->zk[ip];
    }
  }else{
    for(ip=0; ip<np; ip++){
      mtau[2*ip    ] = rho[2*ip    ]*ked1->zk[ip];
      mtau[2*ip + 1] = rho[2*ip + 1]*ked2->zk[ip];
    }
  }

  xc_mgga_evaluate_functional_new
    (p->func_aux[0], max_order, np, rho, sigma, lapl, mtau, mgga);

  /* now we have to combine the results */
  if(out->zk != NULL)
    for(ip=0; ip<np; ip++){
      out->VAR(zk, ip, 0) = mgga->VAR(zk, ip, 0);
    }
#ifndef XC_DONT_COMPILE_VXC
  if(out->vrho != NULL)
    for(ip=0; ip<np; ip++){
#include "maple2c/deorbitalize_1.c"
    }
#ifndef XC_DONT_COMPILE_FXC
  if(out->v2rho2 != NULL)
    for(ip=0; ip<np; ip++){
#include "maple2c/deorbitalize_2.c"
    }
#ifndef XC_DONT_COMPILE_KXC
  if(out->v3rho3 != NULL)
    for(ip=0; ip<np; ip++){
#include "maple2c/deorbitalize_3.c"
    }
#ifndef XC_DONT_COMPILE_LXC
  if(out->v4rho4 != NULL)
    for(ip=0; ip<np; ip++){
#include "maple2c/deorbitalize_4.c"
    }
#endif
#endif
#endif
#endif
  
  xc_deallocate_output_variables(mgga);
  xc_deallocate_output_variables(ked1);
  free(mtau);
  if(p->nspin == XC_POLARIZED){
    xc_deallocate_output_variables(ked2);
    free(mrho); free(msigma); free(mlapl);
  }
}
  
xc_mgga_funcs_variants xc_deorbitalize_func =
  {
   {deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work},
   {deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work, deorbitalized_functional_work},
  };
                                               
