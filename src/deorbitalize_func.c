/*
  Copyright (C) 2006-2007 M.A.L. Marques
                2018-2019 Susi Lehtola
                2019 X. Andrade

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

void xc_evaluate_func(const xc_func_type *func, int max_order,
                      const xc_input_variables *in, xc_output_variables *out)
{
  int ii, check;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  /* turn off orders larger than max_order */
  for(ii=max_order+1; ii <= XC_MAXIMUM_ORDER; ii++)
    orders[ii] = 0;

  /* check if all variables make sense */
  check = xc_input_variables_sanity_check(in, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    fprintf(stderr, "Field %s is not allocated\n", xc_input_variables_name[check]);
    exit(1);
  }

  check = xc_output_variables_sanity_check(out, orders, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    if(check >= 1000)
      fprintf(stderr, "Functional does not provide an implementation of the %d-th derivative\n", check-1000);
    else
      fprintf(stderr, "Field %s is not allocated\n", xc_output_variables_name[check]);
    exit(1);
  }

  xc_output_variables_initialize(out, in->np, func->nspin);

  /* call the specific work routine */
  if(func->info->work != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->work->unpol[max_order] != NULL)
        func->info->work->unpol[max_order](func, in, out);
    }else{
      if(func->info->work->pol[max_order] != NULL)
        func->info->work->pol[max_order](func, in, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, in, out);
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
deorb_work(const xc_func_type *p,
           const xc_input_variables *in, xc_output_variables *out)
{
  double *mtau, *mrho, *msigma, *mlapl;
  xc_input_variables  *in2;
  xc_output_variables *mgga, *ked1, *ked2;
  size_t ip;

  int ii, max_order, flags;
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

  /* prepare buffers for input variables */
  in2 = xc_input_variables_allocate(-1, XC_FAMILY_MGGA, 0, 0);
  in2->np = in->np;

  mtau = (double *) libxc_malloc(in->np*p->dim->tau*sizeof(double));
  if(p->nspin == XC_POLARIZED){
    mrho   = (double *) libxc_malloc(in->np*p->dim->rho*sizeof(double));
    msigma = (double *) libxc_malloc(in->np*p->dim->sigma*sizeof(double));
    mlapl  = (double *) libxc_malloc(in->np*p->dim->lapl*sizeof(double));
  }

  /*
     prepare buffers that will hold the output from the individual functionals
     we have to declare them as XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU
     as the code accesses these values (even if they are zero)
  */
  flags = p->info->flags | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_NEEDS_TAU;
  mgga = xc_output_variables_allocate(in->np, orders, XC_FAMILY_MGGA, flags, p->nspin);
  xc_output_variables_initialize(mgga, in->np, p->nspin);

  ked1 = xc_output_variables_allocate(in->np, orders, XC_FAMILY_MGGA, flags, p->nspin);
  xc_output_variables_initialize(ked1, in->np, p->nspin);
  if(p->nspin == XC_POLARIZED){
    ked2 = xc_output_variables_allocate(in->np, orders, XC_FAMILY_MGGA, flags, p->nspin);
    xc_output_variables_initialize(ked2, in->np, p->nspin);
  }

  /* evaluate the kinetic energy functional */
  if(p->nspin == XC_UNPOLARIZED){
    xc_evaluate_func(p->func_aux[1], max_order, in, ked1);
  }else{
    in2->rho = mrho; in2->sigma = msigma;
    in2->lapl = mlapl;

    for(ip=0; ip<in->np; ip++){
      mrho  [2*ip] = in->rho  [2*ip]; mrho  [2*ip+1] = 0.0;
      msigma[3*ip] = in->sigma[3*ip]; msigma[3*ip+1] = 0.0; msigma[3*ip+2] = 0.0;
      mlapl [2*ip] = in->lapl [2*ip]; mlapl [2*ip+1] = 0.0;
    }
    xc_evaluate_func(p->func_aux[1], max_order, in2, ked1);

    for(ip=0; ip<in->np; ip++){
      mrho  [2*ip] = in->rho  [2*ip + 1];
      msigma[3*ip] = in->sigma[3*ip + 2];
      mlapl [2*ip] = in->lapl [2*ip + 1];
    }
    xc_evaluate_func(p->func_aux[1], max_order, in2, ked2);
  }

  /* now evaluate the mgga functional */
  in2->rho = in->rho; in2->sigma = in->sigma;
  in2->lapl = in->lapl; in2->tau = mtau;
  if(p->nspin == XC_UNPOLARIZED){
    for(ip=0; ip<in->np; ip++){
      mtau[ip] = in->rho[ip]*ked1->zk[ip];
    }
  }else{
    for(ip=0; ip<in->np; ip++){
      mtau[2*ip    ] = in->rho[2*ip    ]*ked1->zk[ip];
      mtau[2*ip + 1] = in->rho[2*ip + 1]*ked2->zk[ip];
    }
  }

  xc_evaluate_func(p->func_aux[0], max_order, in2, mgga);

  /* now we have to combine the results */
  if(out->zk != NULL)
    for(ip=0; ip<in->np; ip++){
      out->VAR(zk, ip, 0) = mgga->VAR(zk, ip, 0);
    }
#ifndef XC_DONT_COMPILE_VXC
  if(out->vrho != NULL)
    for(ip=0; ip<in->np; ip++){
#include "maple2c/deorbitalize_1.c"
    }
#ifndef XC_DONT_COMPILE_FXC
  if(out->v2rho2 != NULL)
    for(ip=0; ip<in->np; ip++){
#include "maple2c/deorbitalize_2.c"
    }
#ifndef XC_DONT_COMPILE_KXC
  if(out->v3rho3 != NULL)
    for(ip=0; ip<in->np; ip++){
#include "maple2c/deorbitalize_3.c"
    }
#ifndef XC_DONT_COMPILE_LXC
  if(out->v4rho4 != NULL)
    for(ip=0; ip<in->np; ip++){
#include "maple2c/deorbitalize_4.c"
    }
#endif
#endif
#endif
#endif

  /* clean up */
  xc_output_variables_deallocate(mgga);
  xc_output_variables_deallocate(ked1);
  libxc_free(in2);
  libxc_free(mtau);
  if(p->nspin == XC_POLARIZED){
    xc_output_variables_deallocate(ked2);
    libxc_free(mrho); libxc_free(msigma); libxc_free(mlapl);
  }
}

xc_functionals_work_variants xc_deorbitalize_func =
  {
   {deorb_work, deorb_work, deorb_work, deorb_work, deorb_work},
   {deorb_work, deorb_work, deorb_work, deorb_work, deorb_work},
  };
