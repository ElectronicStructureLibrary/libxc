/*
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <stdio.h>

#include <xc.h>

#ifdef HAVE_CUDA
#include <cuda.h>
#endif

int main()
{

  int N = 5;
  double rho[N] = {0.1, 0.2, 0.3, 0.4, 0.5};
  double sigma[N] = {0.2, 0.3, 0.4, 0.5, 0.6};
  double exc[N];
  int i, vmajor, vminor, vmicro, func_id = 1;

  xc_version(&vmajor, &vminor, &vmicro);
  printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

  xc_func_type * func = xc_func_alloc();
  
  if(xc_func_init(func, func_id, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", func_id);
    return 1;
  }

  double * rho_cuda;
  cudaMalloc(&rho_cuda, N*sizeof(double));

  double * exc_cuda;
  cudaMalloc(&exc_cuda, N*sizeof(double));

  cudaMemcpy(rho_cuda, rho, N*sizeof(double), cudaMemcpyHostToDevice);
  
  switch(func->info->family)
  {
  case XC_FAMILY_LDA:
    xc_lda_exc(func, 5, rho_cuda, exc_cuda);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga_exc(func, 5, rho, sigma, exc);
    break;
  }

  cudaMemcpy(exc, exc_cuda, N*sizeof(double), cudaMemcpyDeviceToHost);
  
  for(i=0; i<5; i+=1){
    printf("%lf %le\n", rho[i], exc[i]);
  }

  xc_func_end(func);
  xc_func_free(func);
}
