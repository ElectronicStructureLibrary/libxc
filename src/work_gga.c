/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_gga.c
 * @brief This file is to be included in GGA functionals.
 */

/* define auxiliary functions to NULL in case they are not available */
#if defined(XC_DONT_COMPILE_EXC) || maple2c_order < 0 || defined(XC_NO_EXC)
#define work_gga_exc_unpol NULL
#define work_gga_exc_pol NULL
#else
#define ORDER_TXT exc
#define SPIN_TXT  unpol
#include "work_gga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_gga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_VXC) || maple2c_order < 1
#define work_gga_vxc_unpol NULL
#define work_gga_vxc_pol NULL
#else
#define OUT_PARAMS XC_COMMA zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA, )
#endif

#if defined(XC_DONT_COMPILE_FXC) || maple2c_order < 2
#define work_gga_fxc_unpol NULL
#define work_gga_fxc_pol NULL
#else
#define ORDER_TXT fxc
#define SPIN_TXT  unpol
#include "work_gga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_gga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_KXC) || maple2c_order < 3
#define work_gga_kxc_unpol NULL
#define work_gga_kxc_pol NULL
#else
#define ORDER_TXT kxc
#define SPIN_TXT  unpol
#include "work_gga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_gga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_LXC) || maple2c_order < 4
#define work_gga_lxc_unpol NULL
#define work_gga_lxc_pol NULL
#else
  size_t ip;
  double dens;
  double my_rho[2] = {0.0, 0.0};
  double my_sigma[3] = {0.0, 0.0, 0.0};

  for(ip = 0; ip < np; ip++){
    /* Screen low density */
    dens = (p->nspin == XC_POLARIZED) ? rho[0] + rho[1] : rho[0];
    if(dens >= p->dens_threshold) {
      /* sanity check of input parameters */
      my_rho[0] = m_max(p->dens_threshold, rho[0]);
      my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[0]);
      if(p->nspin == XC_POLARIZED){
        double s_ave;

        my_rho[1] = m_max(p->dens_threshold, rho[1]);
        my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[2]);

        my_sigma[1] = sigma[1];
        s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
        /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
        my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
        /* Since |grad n_up - grad n_down|^2 > 0 we also have */
        my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
      }

      if(p->nspin == XC_UNPOLARIZED)
        func_unpol(p, order, my_rho, my_sigma OUT_PARAMS);
      else if(p->nspin == XC_POLARIZED)
        func_pol  (p, order, my_rho, my_sigma OUT_PARAMS);
    }

    /* check for NaNs */
#ifdef XC_DEBUG
    {
      size_t ii;
      const xc_dimensions *dim = &(p->dim);
      int is_OK = 1;

      if(zk != NULL)
        is_OK = is_OK & isfinite(*zk);

      if(vrho != NULL){
        for(ii=0; ii < dim->vrho; ii++)
          is_OK = is_OK && isfinite(vrho[ii]);
        for(ii=0; ii < dim->vsigma; ii++)
          is_OK = is_OK && isfinite(vsigma[ii]);
      }

      if(!is_OK){
        printf("Problem in the evaluation of the functional\n");
        if(p->nspin == XC_UNPOLARIZED){
          printf("./xc-get_data %d 1 %le 0.0 %le 0.0 0.0 0.0 0.0 0.0 0.0\n",
                 p->info->number, *rho, *sigma);
        }else{
          printf("./xc-get_data %d 2 %le %le %le %le %le 0.0 0.0 0.0 0.0\n",
                 p->info->number, rho[0], rho[1], sigma[0], sigma[1], sigma[2]);
        }
      }
    }
#endif

    internal_counters_gga_next(&(p->dim), 0, &rho, &sigma,
                               &zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));
  }   /* for(ip) */

#endif
}

#ifdef HAVE_CUDA

__global__ static void
work_gga_gpu(const XC(func_type) *p, int order, size_t np, const double *rho, const double *sigma,
             double *zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{
  size_t ip = blockIdx.x*blockDim.x + threadIdx.x;
  double dens;
  double my_rho[2] = {0.0, 0.0};
  double my_sigma[3] = {0.0, 0.0, 0.0};

  if(ip >= np) return;

  internal_counters_gga_random(&(p->dim), ip, 0, &rho, &sigma,
                               &zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA &, ));

  /* Density screening */
  dens = (p->nspin == XC_POLARIZED) ? rho[0]+rho[1] : rho[0];
  if(dens >= p->dens_threshold) {
    /* sanity check on input parameters */
    my_rho[0]   = m_max(p->dens_threshold, rho[0]);
    my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[0]);
    if(p->nspin == XC_POLARIZED){
      double s_ave;
      
      my_rho[1]   = m_max(p->dens_threshold, rho[1]);
      my_sigma[2] = m_max(p->sigma_threshold * p->sigma_threshold, sigma[2]);
      
      my_sigma[1] = sigma[1];
      s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
      /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
      my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
      /* Since |grad n_up - grad n_down|^2 > 0 we also have */
      my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
    }
    
    if(p->nspin == XC_UNPOLARIZED)
      func_unpol(p, order, my_rho, my_sigma OUT_PARAMS);
    else if(p->nspin == XC_POLARIZED)
      func_pol  (p, order, my_rho, my_sigma OUT_PARAMS);
  }
}

#endif
