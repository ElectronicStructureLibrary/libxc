/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_mgga.c"
#include "funcs_hyb_mgga.c"

void 
xc_mgga(const xc_func_type *func, int np,
        const double *rho, const double *sigma, const double *lapl, const double *tau,
        double *zk, MGGA_OUT_PARAMS_NO_EXC(double *))
{
  assert(func != NULL);
  const xc_dimensions *dim = &(func->dim);
  
  /* sanity check */
  if(zk != NULL && !(func->info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
            func->info->name);
    exit(1);
  }

  if(vrho != NULL && !(func->info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
            func->info->name);
    exit(1);
  }

  if(v2rho2 != NULL && !(func->info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
            func->info->name);
    exit(1);
  }

  if(v3rho3 != NULL && !(func->info->flags & XC_FLAGS_HAVE_KXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
            func->info->name);
    exit(1);
  }


  /* initialize output to zero */
  if(zk != NULL)
    memset(zk, 0, dim->zk*np*sizeof(double));

  if(vrho != NULL){
    assert(vsigma != NULL);
    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
      assert(vlapl != NULL);
    assert(vtau   != NULL);

    memset(vrho,   0, dim->vrho  *np*sizeof(double));
    memset(vsigma, 0, dim->vsigma*np*sizeof(double));
    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN)
       memset(vlapl,  0, dim->vlapl *np*sizeof(double));
    memset(vtau,   0, dim->vtau  *np*sizeof(double));
  }

  if(v2rho2 != NULL){
    assert(v2rhosigma != NULL);
    assert(v2rhotau   != NULL);
    assert(v2sigma2   != NULL);
    assert(v2sigmatau != NULL);
    assert(v2tau2     != NULL);

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      assert(v2rholapl   != NULL);
      assert(v2sigmalapl != NULL);
      assert(v2lapl2     != NULL);
      assert(v2lapltau   != NULL);
    }
      
    memset(v2rho2,     0, dim->v2rho2     *np*sizeof(double));
    memset(v2rhosigma, 0, dim->v2rhosigma *np*sizeof(double));
    memset(v2rhotau,   0, dim->v2rhotau   *np*sizeof(double));
    memset(v2sigma2,   0, dim->v2sigma2   *np*sizeof(double));
    memset(v2sigmatau, 0, dim->v2sigmatau *np*sizeof(double));
    memset(v2tau2,     0, dim->v2tau2     *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      memset(v2rholapl,   0, dim->v2rholapl  *np*sizeof(double));
      memset(v2sigmalapl, 0, dim->v2sigmalapl*np*sizeof(double));
      memset(v2lapl2,     0, dim->v2lapl2    *np*sizeof(double));
      memset(v2lapltau,   0, dim->v2lapltau  *np*sizeof(double));
    }
  }

  if(v3rho3 != NULL){
    assert(v3rho2sigma   != NULL);
    assert(v3rho2tau     != NULL);
    assert(v3rhosigma2   != NULL);
    assert(v3rhosigmatau != NULL);
    assert(v3rhotau2     != NULL);
    assert(v3sigma3      != NULL);
    assert(v3sigma2tau   != NULL);
    assert(v3sigmatau2   != NULL);
    assert(v3tau3        != NULL);

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){    
      assert(v3rho2lapl     != NULL);
      assert(v3rhosigmalapl != NULL);
      assert(v3rholapl2     != NULL);
      assert(v3rholapltau   != NULL);
      assert(v3sigma2lapl   != NULL);
      assert(v3sigmalapl2   != NULL);
      assert(v3sigmalapltau != NULL);
      assert(v3lapl3        != NULL);
      assert(v3lapl2tau     != NULL);
    }
	
    memset(v3rho3,        0, dim->v3rho3       *np*sizeof(double));
    memset(v3rho2sigma,   0, dim->v3rho2sigma  *np*sizeof(double));
    memset(v3rho2tau,     0, dim->v3rho2tau    *np*sizeof(double));
    memset(v3rhosigma2,   0, dim->v3rhosigma2  *np*sizeof(double));
    memset(v3rhosigmatau, 0, dim->v3rhosigmatau*np*sizeof(double));
    memset(v3rhotau2,     0, dim->v3rhotau2    *np*sizeof(double));
    memset(v3sigma3,      0, dim->v3sigma3     *np*sizeof(double));
    memset(v3sigma2tau,   0, dim->v3sigma2tau  *np*sizeof(double));
    memset(v3sigmatau2,   0, dim->v3sigmatau2  *np*sizeof(double));  
    memset(v3tau3,        0, dim->v3tau3       *np*sizeof(double));

    if(func->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
      memset(v3rho2lapl,     0, dim->v3rho2lapl    *np*sizeof(double));
      memset(v3rhosigmalapl, 0, dim->v3rhosigmalapl*np*sizeof(double));      
      memset(v3rholapl2,     0, dim->v3rholapl2    *np*sizeof(double));
      memset(v3rholapltau,   0, dim->v3rholapltau  *np*sizeof(double));
      memset(v3sigma2lapl,   0, dim->v3sigma2lapl  *np*sizeof(double));
      memset(v3sigmalapl2,   0, dim->v3sigmalapl2  *np*sizeof(double));
      memset(v3sigmalapltau, 0, dim->v3sigmalapltau*np*sizeof(double));
      memset(v3lapl3,        0, dim->v3lapl3       *np*sizeof(double));
      memset(v3lapl2tau,     0, dim->v3lapl2tau    *np*sizeof(double));
    }
  }

  /* call functional */
  if(func->info->mgga != NULL)
    func->info->mgga(func, np, rho, sigma, lapl, tau, zk, MGGA_OUT_PARAMS_NO_EXC());

  /* WARNING: Kxc is not properly mixed */
  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, sigma, lapl, tau,
                zk, vrho, vsigma, vlapl, vtau,
                v2rho2, v2rhosigma, v2rholapl, v2rhotau, 
                v2sigma2, v2sigmalapl, v2sigmatau,
                v2lapl2, v2lapltau,
                v2tau2,
                v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau,
                v3rhosigma2, v3rhosigmalapl, v3rhosigmatau,
                v3rholapl2, v3rholapltau,
                v3rhotau2,
                v3sigma3, v3sigma2lapl, v3sigma2tau,
                v3sigmalapl2, v3sigmalapltau,
                v3sigmatau2,
                v3lapl3, v3lapl2tau,
                v3lapltau2,
                v3tau3
                );

}

/* specializations */
void
xc_mgga_exc(const xc_func_type *p, int np, 
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *zk)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, NULL, NULL, NULL, NULL, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void
xc_mgga_exc_vxc(const xc_func_type *p, int np,
                const double *rho, const double *sigma, const double *lapl, const double *tau,
                double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void
xc_mgga_vxc(const xc_func_type *p, int np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga(p, np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
}

void
xc_mgga_fxc(const xc_func_type *p, int np,
            const double *rho, const double *sigma, const double *lapl, const double *tau,
            double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
            double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
            double *v2lapl2, double *v2lapltau,
            double *v2tau2)
{
  xc_mgga(p, np, rho, sigma, lapl, tau,
          NULL, NULL, NULL, NULL, NULL,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau,
          v2sigma2, v2sigmalapl, v2sigmatau,
          v2lapl2, v2lapltau,
          v2tau2,
          NULL, NULL, NULL, NULL,
          NULL, NULL, NULL,
          NULL, NULL,
          NULL,
          NULL, NULL, NULL,
          NULL, NULL,
          NULL,
          NULL, NULL,
          NULL,
          NULL);
}

void xc_mgga_kxc(const xc_func_type *p, int np,
                 const double *rho, const double *sigma, const double *lapl, const double *tau,
                 double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                 double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                 double *v3rholapl2, double *v3rholapltau,
                 double *v3rhotau2,
                 double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
                 double *v3sigmalapl2, double *v3sigmalapltau,
                 double *v3sigmatau2,
                 double *v3lapl3, double *v3lapl2tau,
                 double *v3lapltau2,
                 double *v3tau3)
{
  xc_mgga(p, np, rho, sigma, lapl, tau,
          NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL,
          NULL, NULL, NULL,
          NULL, NULL,
          NULL,
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau,
          v3rhosigma2, v3rhosigmalapl, v3rhosigmatau,
          v3rholapl2, v3rholapltau,
          v3rhotau2,
          v3sigma3, v3sigma2lapl, v3sigma2tau,
          v3sigmalapl2, v3sigmalapltau,
          v3sigmatau2,
          v3lapl3, v3lapl2tau,
          v3lapltau2,
          v3tau3);
}
