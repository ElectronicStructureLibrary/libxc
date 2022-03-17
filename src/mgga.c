/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_mgga.c"

void xc_evaluate_mgga(const xc_func_type *func, int order, size_t np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                xc_output_variables *out)
{
  int ii, check;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  /* turn off orders smaller than max_order */
  for(ii=max_order+1; ii <= XC_MAXIMUM_ORDER; ii++)
    orders[ii] = 0;

  /* check if all variables make sense */
  check = xc_output_variables_sanity_check(out, orders, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    if(check >= 1000)
      fprintf(stderr, "Functional does not provide an implementation of the %d-th derivative\n", check-1000);
    else
      fprintf(stderr, "Field %s is not allocated\n", xc_output_variables_name[check]);
    exit(1);
  }
  
  xc_output_variables_initialize(out, np, func->nspin);

  /* call the mGGA routines */
  if(func->info->mgga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->mgga->unpol[order] != NULL)
        func->info->mgga->unpol[order](func, np, rho, sigma, lapl, tau, out);
    }else{
      if(func->info->mgga->pol[order] != NULL)
        func->info->mgga->pol[order](func, np, rho, sigma, lapl, tau, out);
    }
  }
    
  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, lapl, tau, NULL, out);
}

/* old API */
#define SET_ORDER_0                          \
  out.zk = zk;

#define SET_ORDER_1                          \
  out.vrho = vrho;                           \
  out.vsigma = vsigma;                       \
  out.vlapl = vlapl;                         \
  out.vtau = vtau;

#define SET_ORDER_2                          \
  out.v2rho2 = v2rho2;                       \
  out.v2rhosigma = v2rhosigma;               \
  out.v2rholapl = v2rholapl;                 \
  out.v2rhotau = v2rhotau;                   \
  out.v2sigma2 = v2sigma2;                   \
  out.v2sigmalapl = v2sigmalapl;             \
  out.v2sigmatau = v2sigmatau;               \
  out.v2lapl2 = v2lapl2;                     \
  out.v2lapltau = v2lapltau;                 \
  out.v2tau2 = v2tau2;

#define SET_ORDER_3                          \
  out.v3rho3 = v3rho3;                       \
  out.v3rho2sigma = v3rho2sigma;             \
  out.v3rho2lapl = v3rho2lapl;               \
  out.v3rho2tau = v3rho2tau;                 \
  out.v3rhosigma2 = v3rhosigma2;             \
  out.v3rhosigmalapl = v3rhosigmalapl;       \
  out.v3rhosigmatau = v3rhosigmatau;         \
  out.v3rholapl2 = v3rholapl2;               \
  out.v3rholapltau = v3rholapltau;           \
  out.v3rhotau2 = v3rhotau2;                 \
  out.v3sigma3 = v3sigma3;                   \
  out.v3sigma2lapl = v3sigma2lapl;           \
  out.v3sigma2tau = v3sigma2tau;             \
  out.v3sigmalapl2 = v3sigmalapl2;           \
  out.v3sigmalapltau = v3sigmalapltau;       \
  out.v3sigmatau2 = v3sigmatau2;             \
  out.v3lapl3 = v3lapl3;                     \
  out.v3lapl2tau = v3lapl2tau;               \
  out.v3lapltau2 = v3lapltau2;               \
  out.v3tau3 = v3tau3;

#define SET_ORDER_4                          \
  out.v4rho4 = v4rho4;                       \
  out.v4rho3sigma = v4rho3sigma;             \
  out.v4rho3lapl = v4rho3lapl;               \
  out.v4rho3tau = v4rho3tau;                 \
  out.v4rho2sigma2 = v4rho2sigma2;           \
  out.v4rho2sigmalapl = v4rho2sigmalapl;     \
  out.v4rho2sigmatau = v4rho2sigmatau;       \
  out.v4rho2lapl2 = v4rho2lapl2;             \
  out.v4rho2lapltau = v4rho2lapltau;         \
  out.v4rho2tau2 = v4rho2tau2;               \
  out.v4rhosigma3 = v4rhosigma3;             \
  out.v4rhosigma2lapl = v4rhosigma2lapl;     \
  out.v4rhosigma2tau = v4rhosigma2tau;       \
  out.v4rhosigmalapl2 = v4rhosigmalapl2;     \
  out.v4rhosigmalapltau = v4rhosigmalapltau; \
  out.v4rhosigmatau2 = v4rhosigmatau2;       \
  out.v4rholapl3 = v4rholapl3;               \
  out.v4rholapl2tau = v4rholapl2tau;         \
  out.v4rholapltau2 = v4rholapltau2;         \
  out.v4rhotau3 = v4rhotau3;                 \
  out.v4sigma4 = v4sigma4;                   \
  out.v4sigma3lapl = v4sigma3lapl;           \
  out.v4sigma3tau = v4sigma3tau;             \
  out.v4sigma2lapl2 = v4sigma2lapl2;         \
  out.v4sigma2lapltau = v4sigma2lapltau;     \
  out.v4sigma2tau2 = v4sigma2tau2;           \
  out.v4sigmalapl3 = v4sigmalapl3;           \
  out.v4sigmalapl2tau = v4sigmalapl2tau;     \
  out.v4sigmalapltau2 = v4sigmalapltau2;     \
  out.v4sigmatau3 = v4sigmatau3;             \
  out.v4lapl4 = v4lapl4;                     \
  out.v4lapl3tau = v4lapl3tau;               \
  out.v4lapl2tau2 = v4lapl2tau2;             \
  out.v4lapltau3 = v4lapltau3;               \
  out.v4tau4 = v4tau4;

void
xc_mgga(const xc_func_type *p, size_t np,
        const double *rho, const double *sigma, const double *lapl, const double *tau,
        double *zk,
        double *vrho, double *vsigma, double *vlapl, double *vtau,
        double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau, double *v2sigma2,
        double *v2sigmalapl, double *v2sigmatau, double *v2lapl2, double *v2lapltau, double *v2tau2,
        double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau, double *v3rhosigma2,
        double *v3rhosigmalapl, double *v3rhosigmatau, double *v3rholapl2, double *v3rholapltau,
        double *v3rhotau2, double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
        double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2, double *v3lapl3,
        double *v3lapl2tau, double *v3lapltau2, double *v3tau3,
        double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
        double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
        double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
        double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
        double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
        double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
        double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
        double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
        double *v4lapl2tau2, double *v4lapltau3, double *v4tau4
        )
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));

  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;
  SET_ORDER_4;

  /* call new API */
  xc_evaluate_mgga(p, order, np, rho, sigma, lapl, tau, &out);
}


/* specializations */
void
xc_mgga_exc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_0;

  xc_evaluate_mgga(p, 0, np, rho, sigma, lapl, tau, &out);
}

void
xc_mgga_exc_vxc(const xc_func_type *p, size_t np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_0;
  SET_ORDER_1;

  xc_evaluate_mgga(p, 1, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;

  xc_evaluate_mgga(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_1;
  SET_ORDER_2;

  xc_evaluate_mgga(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_0;
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  xc_evaluate_mgga(p, 3, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2,
                         double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                         double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                         double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                         double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                         double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                         double *v3tau3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_1;
  SET_ORDER_2;
  SET_ORDER_3;

  xc_evaluate_mgga(p, 3, np, rho, sigma, lapl, tau, &out);
}


void
xc_mgga_vxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_1;

  xc_evaluate_mgga(p, 1, np, rho, sigma, lapl, tau, &out);
}

void
xc_mgga_fxc(const xc_func_type *p, size_t np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
            double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
            double *v2lapltau, double *v2tau2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_2;

  xc_evaluate_mgga(p, 2, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_kxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                 double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                 double *v3rholapl2, double *v3rholapltau,  double *v3rhotau2,
                 double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
                 double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2,
                 double *v3lapl3, double *v3lapl2tau, double *v3lapltau2, double *v3tau3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_3;

  xc_evaluate_mgga(p, 3, np, rho, sigma, lapl, tau, &out);
}

void xc_mgga_lxc(const xc_func_type *p, size_t np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
                 double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
                 double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
                 double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
                 double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
                 double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
                 double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
                 double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
                 double *v4lapl2tau2, double *v4lapltau3, double *v4tau4)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  SET_ORDER_4;

  xc_evaluate_mgga(p, 4, np, rho, sigma, lapl, tau, &out);
}
