/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"

#ifdef HAVE_FORTRAN

#include "xc.h"
#include "string_f.h"

/* version */
void FC_FUNC(xc_f90_version, XC_F90_VERSION)
     (int *major, int *minor, int *micro)
{
  xc_version(major, minor, micro);
}

void FC_FUNC(xc_f90_version_string, XC_F90_VERSIN_STRING)
     (STR_F_TYPE version_string STR_ARG1)
{
  const char *version;

  version = xc_version_string();
  TO_F_STR1(version, version_string);
}

/* info */
CC_FORTRAN_INT FC_FUNC(xc_f90_info_number, XC_F90_INFO_NUMBER)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->number;
}


CC_FORTRAN_INT FC_FUNC(xc_f90_info_kind, XC_F90_INFO_KIND)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->kind;
}


void FC_FUNC(xc_f90_info_name, XC_F90_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((xc_func_info_type *)(*info))->name, s);
}


CC_FORTRAN_INT  FC_FUNC(xc_f90_info_family, XC_F90_INFO_FAMILY)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->family;
}


CC_FORTRAN_INT  FC_FUNC(xc_f90_info_flags, XC_F90_INFO_FLAGS)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->flags;
}


void FC_FUNC(xc_f90_info_refs, XC_F90_INFO_REFS)
     (void **info, CC_FORTRAN_INT *number, STR_F_TYPE ref_f STR_ARG1)
{
  xc_func_info_type *func_p = (xc_func_info_type *)(*info);

  assert(*number >=0 && *number < 5);

  if(func_p->refs[*number] == NULL){
    *number = -1;
    return;
  }

  TO_F_STR1(func_p->refs[*number]->ref, ref_f);

  (*number)++;
  fflush(stdout);
}


void FC_FUNC(xc_f90_functional_get_name, XC_F90_FUNCTIONAL_GET_NAME)
     (CC_FORTRAN_INT *func_number, STR_F_TYPE func_string STR_ARG1)
{
  char *name;

  name = xc_functional_get_name(*func_number);
  if ( name == NULL ) name = strndup("unknown", 256);

  TO_F_STR1(name, func_string);
  libxc_free(name);
}


CC_FORTRAN_INT  FC_FUNC(xc_f90_functional_get_number, XC_F90_FUNCTIONAL_GET_NUMBER)
     (STR_F_TYPE func_string STR_ARG1)
{
  char *name;
  int ret;

  TO_C_STR1(func_string, name);
  
  ret = xc_functional_get_number(name);
  libxc_free(name);

  return (CC_FORTRAN_INT) ret;
}


/* functionals */
CC_FORTRAN_INT  FC_FUNC(xc_f90_family_from_id, XC_F90_FAMILY_FROM_ID)
  (CC_FORTRAN_INT  *functional)
{
  return (CC_FORTRAN_INT) xc_family_from_id((int) (*functional), NULL, NULL);
}

CC_FORTRAN_INT FC_FUNC(xc_f90_number_of_functionals, XC_F90_NUMBER_OF_FUNCTIONALS)
  ()
{
  return (CC_FORTRAN_INT) xc_number_of_functionals();
}

CC_FORTRAN_INT FC_FUNC(xc_f90_maximum_name_length, XC_F90_MAXIMUM_LENGTH_NAME)
  ()
{
  return (CC_FORTRAN_INT) xc_maximum_name_length();
}

void FC_FUNC(xc_f90_available_functional_numbers, XC_F90_AVAILABLE_FUNCTIONAL_NUMBERS)
  (CC_FORTRAN_INT *list)
{
  xc_available_functional_numbers(list);
}


/* Standard initialization */
void FC_FUNC(xc_f90_func_init, XC_F90_FUNC_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_func_type *func_p;
  
  func_p = (xc_func_type *)libxc_malloc(sizeof(xc_func_type));
  xc_func_init(func_p, (int) (*functional), (int) (*nspin));

  *p    = (void *) func_p;
  *info = (void *)(func_p->info);
}

void FC_FUNC(xc_f90_func_end, XC_F90_FUNC_END)
     (void **p)
{
  xc_func_end((xc_func_type *)(*p));
  libxc_free(*p);
  *p = NULL;
}

void FC_FUNC(xc_f90_func_set_dens_threshold, XC_F90_FUNC_SET_DENS_THRESHOLD)
     (void **p, double *dens_threshold)
{
  xc_func_set_dens_threshold((xc_func_type *)(*p), *dens_threshold);
}

void FC_FUNC(xc_f90_func_set_ext_params, XC_F90_FUNC_SET_EXT_PARAMS)
     (void **p, double *ext_params)
{
  xc_func_set_ext_params((xc_func_type *)(*p), ext_params);
}


/* LDAs */

void FC_FUNC(xc_f90_lda, XC_F90_LDA)
     (void **p, CC_FORTRAN_INT *np, double *rho, 
      double *zk, double *vrho, double *v2rho2, double *v3rho3, double *v4rho4)
{
  xc_lda((xc_func_type *)(*p), *np, rho, zk, vrho, v2rho2, v3rho3, v4rho4);
}

void FC_FUNC(xc_f90_lda_exc, XC_F90_LDA_EXC)
     (void **p, CC_FORTRAN_INT *np, double *rho,
      double *zk)
{
  xc_lda((xc_func_type *)(*p), *np, rho, zk, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_lda_exc_vxc, XC_F90_LDA_EXC_VXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, 
      double *zk, double *vrho)
{
  xc_lda((xc_func_type *)(*p), *np, rho, zk, vrho, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_lda_vxc, XC_F90_LDA_VXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, 
      double *vrho)
{
  xc_lda((xc_func_type *)(*p), *np, rho, NULL, vrho, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_lda_vxc_fxc, XC_F90_LDA_VXC_FXC)
     (void **p, CC_FORTRAN_INT *np, double *rho,
      double *vrho, double *v2rho2)
{
  xc_lda((xc_func_type *)(*p), *np, rho, NULL, vrho, v2rho2, NULL, NULL);
}

void FC_FUNC(xc_f90_lda_fxc, XC_F90_LDA_FXC)
     (void **p, CC_FORTRAN_INT *np, double *rho,
      double *v2rho2)
{
  xc_lda((xc_func_type *)(*p), *np, rho, NULL, NULL, v2rho2, NULL, NULL);
}

void FC_FUNC(xc_f90_lda_kxc, XC_F90_LDA_KXC)
     (void **p, CC_FORTRAN_INT *np, double *rho,
      double *v3rho3)
{
  xc_lda((xc_func_type *)(*p), *np, rho, NULL, NULL, NULL, v3rho3, NULL);
}

void FC_FUNC(xc_f90_lda_lxc, XC_F90_LDA_LXC)
     (void **p, CC_FORTRAN_INT *np, double *rho,
      double *v4rho4)
{
  xc_lda((xc_func_type *)(*p), *np, rho, NULL, NULL, NULL, NULL, v4rho4);
}

/* GGAs */

void FC_FUNC(xc_f90_gga, XC_F90_GGA)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *zk, GGA_OUT_PARAMS_NO_EXC(double *))
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, zk, GGA_OUT_PARAMS_NO_EXC(XC_NOARG));
}

void FC_FUNC(xc_f90_gga_exc, XC_F90_GGA_EXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *zk)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_exc_vxc, XC_F90_GGA_EXC_VXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *zk, double *vrho, double *vsigma)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_vxc, XC_F90_GGA_VXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *vrho, double *vsigma)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, NULL, vrho, vsigma, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_vxc_fxc, XC_F90_GGA_VXC_FXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *vrho, double *vsigma,
      double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, NULL, vrho, vsigma, 
         v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_fxc, XC_F90_GGA_FXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, NULL, NULL, NULL, 
         v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_kxc, XC_F90_GGA_KXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, NULL, NULL, NULL, 
         NULL, NULL, NULL, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,
         NULL, NULL, NULL, NULL, NULL);
}

void FC_FUNC(xc_f90_gga_lxc, XC_F90_GGA_LXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, 
      double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4)
{
  xc_gga((xc_func_type *)(*p), *np, rho, sigma, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4);
}

void FC_FUNC(xc_f90_gga_ak13_get_asymptotic, XC_F90_GGA_AK13_GET_ASYMPTOTIC)
  (double *homo, double *asymp)
{
  *asymp = xc_gga_ak13_get_asymptotic(*homo);
}

void FC_FUNC(xc_f90_hyb_exx_coef, XC_F90_HYB_EXX_COEF)
   (void **p, double *coef)
{
  *coef = xc_hyb_exx_coef((xc_func_type *)(*p));
}

void FC_FUNC(xc_f90_hyb_cam_coef, XC_F90_HYB_CAM_COEF)
  (void **p, double *omega, double *alpha, double *beta)
{
  xc_hyb_cam_coef((xc_func_type *)(*p), omega, alpha, beta);
}

void FC_FUNC(xc_f90_nlc_coef, XC_F90_NLC_COEF)
  (void **p, double *nlc_b, double *nlc_c)
{
  xc_nlc_coef((xc_func_type *)(*p), nlc_b, nlc_c);
}

/* meta-GGAs */

void FC_FUNC(xc_f90_mgga, XC_F90_MGGA)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
      double *zk, MGGA_OUT_PARAMS_NO_EXC(double *)
      )
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, 
          zk, MGGA_OUT_PARAMS_NO_EXC(XC_NOARG));
}

void FC_FUNC(xc_f90_mgga_exc, XC_F90_MGGA_EXC)
     (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau, 
      double *zk)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, zk, NULL, NULL, NULL, NULL, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void FC_FUNC(xc_f90_mgga_exc_vxc, XC_F90_MGGA_EXC_VXC)
      (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
       double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void FC_FUNC(xc_f90_mgga_vxc, XC_F90_MGGA_VXC)
      (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
       double *vrho, double *vsigma, double *vlapl, double *vtau)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau, 
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void FC_FUNC(xc_f90_mgga_vxc_fxc, XC_F90_MGGA_VXC_FXC)
      (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
       double *vrho, double *vsigma, double *vlapl, double *vtau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau, 
       double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2lapl2, double *v2lapltau,
       double *v2tau2)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, NULL, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
           );
}

void FC_FUNC(xc_f90_mgga_fxc, XC_F90_MGGA_FXC)
      (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau, 
       double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2lapl2, double *v2lapltau,
       double *v2tau2)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}
      
      
void FC_FUNC(xc_f90_mgga_kxc, XC_F90_MGGA_KXC)
  (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
   double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
   double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
   double *v3rholapl2, double *v3rholapltau,  double *v3rhotau2,
   double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
   double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2,
   double *v3lapl3, double *v3lapl2tau, double *v3lapltau2, double *v3tau3)
{
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,
          v3lapltau2, v3tau3,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
          NULL, NULL, NULL, NULL, NULL
          );
}

void FC_FUNC(xc_f90_mgga_lxc, XC_F90_MGGA_LXC)
  (void **p, CC_FORTRAN_INT *np, double *rho, double *sigma, double *lapl, double *tau,
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
  xc_mgga((xc_func_type *)(*p), *np, rho, sigma, lapl, tau, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2,
            v4rho2sigmalapl, v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau,
            v4rho2tau2, v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau,
            v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmatau2,
            v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3,
            v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma2lapl2,
            v4sigma2lapltau, v4sigma2tau2, v4sigmalapl3, v4sigmalapl2tau,
            v4sigmalapltau2, v4sigmatau3, v4lapl4, v4lapl3tau,
            v4lapl2tau2, v4lapltau3, v4tau4
          );
}

#endif
