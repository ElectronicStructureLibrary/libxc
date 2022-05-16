/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_lda_inc.c
 * @brief This file is to be included in LDA functionals.
 */

#ifdef XC_DEBUG
#define __USE_GNU 1
#include <fenv.h>
#endif

/* macro to simplify accessing the variables */
#define VAR(var, ip, index)        var[ip*p->dim->var + index]
#define WORK_LDA_(order, spin)     work_lda_ ## order ## _ ## spin
#define WORK_LDA_IP_(order, spin)  work_lda_ip_ ## order ## _ ## spin
#define FUNC_(order, spin)         func_     ## order ## _ ## spin

/* we need double escaping of the preprocessor macros */
#define WORK_LDA(order, spin)     WORK_LDA_(order, spin)
#define WORK_LDA_IP(order, spin)  WORK_LDA_IP_(order, spin)
#define FUNC(order, spin)         FUNC_(order, spin)

#ifdef HAVE_CUDA
__global__ static void
WORK_LDA_IP(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, const xc_input_variables *in, xc_output_variables *out)
#else
inline static void
WORK_LDA_IP(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, size_t ip, const xc_input_variables *in, xc_output_variables *out)
#endif
{
  double dens; /* the total density */
  double my_rho[2] = {0.0, 0.0};

#ifdef HAVE_CUDA
  size_t ip = blockIdx.x*blockDim.x + threadIdx.x;
#else
#ifdef XC_DEBUG
  /* This throws an exception when floating point errors are encountered */
  /* feenableexcept(FE_DIVBYZERO | FE_INVALID); */
#endif
#endif

  /* this check is required for the GPU kernel */
  if(ip >= in->np)
    return;
  
  /* screen small densities */
  dens = (p->nspin == XC_POLARIZED) ?
    in->VAR(rho, ip, 0) + in->VAR(rho, ip, 1) :
    in->VAR(rho, ip, 0);
  if(dens < p->dens_threshold)
    return;
    
  /* sanity check of input parameters */
  my_rho[0] = m_max(p->dens_threshold, in->VAR(rho, ip, 0));
  if(p->nspin == XC_POLARIZED){
    my_rho[1] = m_max(p->dens_threshold, in->VAR(rho, ip, 1));
  }

  /* evaluate the functional */
  FUNC(ORDER_TXT, SPIN_TXT)(p, ip, my_rho, out);

#ifdef XC_DEBUG
  /* check for NaNs in the output */
  const xc_dimensions *dim = p->dim;
  int ii, is_OK = 1;

  if(out->zk != NULL)
    is_OK = is_OK & isfinite(out->VAR(zk, ip, 0));

  if(out->vrho != NULL){
    for(ii=0; ii < dim->vrho; ii++)
      is_OK = is_OK && isfinite(out->VAR(vrho, ip, ii));
  }
  
  if(!is_OK){
    printf("Problem in the evaluation of the functional\n");
    if(p->nspin == XC_UNPOLARIZED){
      printf("./xc-get_data %d 1 %le 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
             p->info->number, in->VAR(rho, ip, 0));
    }else{
      printf("./xc-get_data %d 2 %le %le 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
             p->info->number, in->VAR(rho, ip, 0), in->VAR(rho, ip, 1));
    }
  }
#endif
}

static void
WORK_LDA(ORDER_TXT, SPIN_TXT)
(const XC(func_type) *p, const xc_input_variables *in, xc_output_variables *out)
{
#ifdef HAVE_CUDA
  /* make a copy of 'p' and 'out' since they might be in host-only memory */
  XC(func_type) *pcuda = (XC(func_type) *) libxc_malloc(sizeof(XC(func_type)));
  xc_input_variables *incuda = (xc_input_variables *) libxc_malloc(sizeof(xc_input_variables));
  xc_output_variables *outcuda = (xc_output_variables *) libxc_malloc(sizeof(xc_output_variables));
  
  cudaMemcpy(pcuda, p, sizeof(XC(func_type)), cudaMemcpyHostToDevice);
  cudaMemcpy(incuda, in, sizeof(xc_input_variables), cudaMemcpyHostToDevice);
  cudaMemcpy(outcuda, out, sizeof(xc_output_variables), cudaMemcpyHostToDevice);

  /* determine number of blocks required */
  size_t nblocks = in->np/CUDA_BLOCK_SIZE;
  if(in->np != nblocks*CUDA_BLOCK_SIZE) nblocks++;

  /* execute kernel */
  WORK_LDA_IP(ORDER_TXT, SPIN_TXT)<<<nblocks, CUDA_BLOCK_SIZE>>>
    (pcuda, incuda, outcuda);

  /* clean up memory */
  libxc_free(pcuda);
  libxc_free(incuda);
  libxc_free(outcuda);
#else
  size_t ip;
  
  /* simply loop over points */
  for(ip=0; ip<in->np; ip++){
    WORK_LDA_IP(ORDER_TXT, SPIN_TXT)(p, ip, in, out);
  }
#endif
}
