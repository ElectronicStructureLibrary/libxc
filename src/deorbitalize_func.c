/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

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
deorb_work(const xc_func_type *p, size_t np,
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
  mgga = xc_output_variables_allocate
    (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);
  xc_output_variables_initialize(mgga, np, p->nspin);
  
  ked1 = xc_output_variables_allocate
    (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);
  xc_output_variables_initialize(ked1, np, p->nspin);

  mtau   = (double *) libxc_malloc(np*p->dim->tau*sizeof(double));
  if(p->nspin == XC_POLARIZED){
    mrho   = (double *) libxc_malloc(np*p->dim->rho*sizeof(double));
    msigma = (double *) libxc_malloc(np*p->dim->sigma*sizeof(double));
    mlapl  = (double *) libxc_malloc(np*p->dim->lapl*sizeof(double));

    ked2 = xc_output_variables_allocate
      (np, orders, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU, p->nspin);
    xc_output_variables_initialize(ked2, np, p->nspin);
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
  
  xc_output_variables_deallocate(mgga);
  xc_output_variables_deallocate(ked1);
  free(mtau);
  if(p->nspin == XC_POLARIZED){
    xc_output_variables_deallocate(ked2);
    free(mrho); free(msigma); free(mlapl);
  }
}
  
xc_mgga_funcs_variants xc_deorbitalize_func =
  {
   {deorb_work, deorb_work, deorb_work, deorb_work, deorb_work},
   {deorb_work, deorb_work, deorb_work, deorb_work, deorb_work},
  };
                                               
