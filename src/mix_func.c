/*
  Copyright (C) 2006-2021 M.A.L. Marques
                2018-2021 Susi Lehtola
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
}

#ifdef HAVE_CUDA
__global__ static void add_to_mix_gpu(size_t np, double * dst, double coeff, const double *src){
  size_t ip = blockIdx.x * blockDim.x + threadIdx.x;
  if(ip < np) dst[ip] += coeff*src[ip];
}
#endif

#define is_hgga(id)   ((id) == XC_FAMILY_HGGA)
#define is_mgga(id)   ((id) == XC_FAMILY_MGGA || is_hgga(id))
#define is_gga(id)    ((id) == XC_FAMILY_GGA  || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA  ||  is_gga(id))

/*
   Sanity check: have we claimed the highest possible derivatives?
   First, check for the lowest common derivative (also need to make
   sure the derivatives have been compiled in!)
*/
void mix_func_sanity_check(const xc_func_type *func)
{
  int ii;
  const xc_func_type *aux;

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

  /* Same for tau */
  int need_tau = 0;
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    if(aux->info->flags & XC_FLAGS_NEEDS_TAU)
      need_tau = XC_FLAGS_NEEDS_TAU;
  }
  assert((func->info->flags & XC_FLAGS_NEEDS_TAU) == need_tau);

  /* Check compatibility of the individual components */
  for(ii=0; ii<func->n_func_aux; ii++){
    aux = func->func_aux[ii];
    /* Sanity check: if component is GGA or meta-GGA, mix functional
       must also be GGA or meta-GGA */
    if(is_gga(aux->info->family))
      assert(is_gga(func->info->family));
    if(is_mgga(aux->info->family) && !is_mgga(func->info->family))
      assert(is_mgga(func->info->family));
    if(is_hgga(aux->info->family) && !is_hgga(func->info->family))
      assert(is_hgga(func->info->family));

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
}

void
xc_mix_func(const xc_func_type *func, const xc_input_variables *in, xc_output_variables *out)
{
  const xc_func_type *aux;
  xc_output_variables *xout;
  size_t ip;

  int ifunc, ii, max_order, check;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  const xc_dimensions *dim = func->dim;

  max_order = -1;
  for(ii=0; ii <= XC_MAXIMUM_ORDER; ii++){
    if(orders[ii]) max_order = ii;
  }

  mix_func_sanity_check(func);

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

  xout = xc_output_variables_allocate
    (in->np, orders, func->info->family, func->info->flags, func->nspin);

  /* Proceed by computing the mix */
  for(ifunc=0; ifunc<func->n_func_aux; ifunc++){
    aux = func->func_aux[ifunc];

    /* have to clean explicitly the buffers here */
    xc_output_variables_initialize(xout, in->np, func->nspin);

    /* Evaluate the functional */
    xc_evaluate_func(aux, max_order, in, xout);

    /* Do the mixing */
    for(ii=0; ii<XC_TOTAL_NUMBER_OUTPUT_VARIABLES; ii++){
      if(out->fields[ii] == NULL)
        continue;
      /* this could be replaced by a daxpy BLAS call */
#ifndef HAVE_CUDA
      for(ip=0; ip<in->np*dim->fields[ii+5]; ip++)
        out->fields[ii][ip] += func->mix_coef[ifunc]*xout->fields[ii][ip];
#else
      size_t nblocks = in->np/CUDA_BLOCK_SIZE;
      if(in->np != nblocks*CUDA_BLOCK_SIZE) nblocks++;
      add_to_mix_gpu<<<nblocks, CUDA_BLOCK_SIZE>>>
        (in->np*dim->fields[ii+5], out->fields[ii], func->mix_coef[ifunc], xout->fields[ii]);
#endif
    }
  }
  xc_output_variables_deallocate(xout);
}

int
xc_num_aux_funcs(const xc_func_type *p) {
  assert(p != NULL);
  return p->n_func_aux;
}

void
xc_aux_func_ids(const xc_func_type *p, int *ids) {
  int i;
  for(i=0; i<p->n_func_aux;i++)
    ids[i] = p->func_aux[i]->info->number;
}

void
xc_aux_func_weights(const xc_func_type *p, double *weights) {
  int i;
  for(i=0; i<p->n_func_aux;i++)
    weights[i] = p->mix_coef[i];
}
