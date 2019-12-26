/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"


/* this function converts the spin-density into total density and
	 relative magnetization */
/* inline */ GPU_FUNCTION void
xc_rho2dzeta(int nspin, const double *rho, double *d, double *zeta)
{
  if(nspin==XC_UNPOLARIZED){
    *d    = max(rho[0], 0.0);
    *zeta = 0.0;
  }else{
    *d = rho[0] + rho[1];
    if(*d > 0.0){
      *zeta = (rho[0] - rho[1])/(*d);
      *zeta = min(*zeta,  1.0);
      *zeta = max(*zeta, -1.0);
    }else{
      *d    = 0.0;
      *zeta = 0.0;
    }
  }
}

const char *get_kind(const xc_func_type *func) {
  switch(func->info->kind) {
   case(XC_EXCHANGE):
      return "XC_EXCHANGE";

    case(XC_CORRELATION):
      return "XC_CORRELATION";

    case(XC_EXCHANGE_CORRELATION):
      return "XC_EXCHANGE_CORRELATION";

    case(XC_KINETIC):
      return "XC_KINETIC";

    default:
      printf("Internal error in get_kind.\n");
      return "";
  }
}

const char *get_family(const xc_func_type *func) {
  switch(func->info->family) {
    case(XC_FAMILY_UNKNOWN):
      return "XC_FAMILY_UNKNOWN";

    case(XC_FAMILY_LDA):
      return "XC_FAMILY_LDA";

    case(XC_FAMILY_GGA):
      return "XC_FAMILY_GGA";

    case(XC_FAMILY_MGGA):
      return "XC_FAMILY_MGGA";

    case(XC_FAMILY_LCA):
      return "XC_FAMILY_LCA";

    case(XC_FAMILY_OEP):
      return "XC_FAMILY_OEP";

    case(XC_FAMILY_HYB_LDA):
      return "XC_FAMILY_HYB_LDA";

    case(XC_FAMILY_HYB_GGA):
      return "XC_FAMILY_HYB_GGA";

    case(XC_FAMILY_HYB_MGGA):
      return "XC_FAMILY_HYB_MGGA";

    default:
      printf("Internal error in get_family.\n");
      return "";
  }
}

/* this function checks if it should use the default or
   the user assigned value for an external parameter */
double
get_ext_param(const func_params_type *params, const double *values, int index)
{
  /* 
     If libxc finds a file in the current directory name
     "libxc.params", it will try to read the parameters for the
     current functional from it. This file should contain one
     parameter per line. E.g., for the x_pbe functional:

       ------------------ <start libxc.params>
       0.8040              # _kappa
       0.2195149727645171  # _mu (PBE)
       ------------------ <end libxc.params>

     Note that this only works for functionals whose parameters can be
     set by set_ext_params.
  */

  /* Commented as considered dangerous ;)
  FILE *par_in;
  int ii, nn;
  double dd;

  if((par_in = fopen("libxc.params","rb"))){
    for(ii=0; ii<index; ii++)
      fscanf(par_in, "%*[^\n]\n", NULL);

    nn = fscanf(par_in, "%lf", &dd);
    fclose(par_in);

    if(nn == 1)
      return dd;
  }
  */

  if(values == NULL || values[index] == XC_EXT_PARAMS_DEFAULT)
    return params[index].value; /* return default value */
  else
    return values[index]; /* return user assigned value */
}

/* these functional handle the internal counters
   used to move along the input and output arrays.
   We have to pay particular attention to the spin,
   of course. */
void
internal_counters_set_lda(int nspin, xc_dimensions *dim)
{
  dim->rho = dim->vrho = nspin;
  dim->zk  = 1;
  if(nspin == XC_UNPOLARIZED){
    dim->v2rho2 = dim->v3rho3 = dim->v4rho4 = 1;
  }else{
    dim->v2rho2 = 3;
    dim->v3rho3 = 4;
    dim->v4rho4 = 5;
  }
}

void
internal_counters_set_gga(int nspin, xc_dimensions *dim)
{
  internal_counters_set_lda(nspin, dim);

  if(nspin == XC_UNPOLARIZED){
    dim->sigma  = dim->vsigma = 1;
    dim->v2rhosigma  = dim->v2sigma2 = 1;
    dim->v3rho2sigma = dim->v3rhosigma2 = dim->v3sigma3 = 1;
    dim->v4rho3sigma = dim->v4rho2sigma2 = dim->v4rhosigma3 = dim->v4sigma4 = 1;
    
  }else{
    dim->sigma = 3;

    dim->vsigma = 3;

    dim->v2rhosigma = 2*3;
    dim->v2sigma2 = 6;
    
    dim->v3rho2sigma = 3*3;
    dim->v3rhosigma2 = 2*6;
    dim->v3sigma3    = 10;

    dim->v4rho3sigma  = 4*3;
    dim->v4rho2sigma2 = 3*6;
    dim->v4rhosigma3  = 2*10;
    dim->v4sigma4     = 15;
  }
}

void
internal_counters_set_mgga(int nspin, xc_dimensions *dim)
{
  internal_counters_set_gga(nspin, dim);

  dim->lapl = dim->vlapl = nspin;
  dim->tau  = dim->vtau  = nspin;
  if(nspin == XC_UNPOLARIZED){
    dim->v2lapl2 = dim->v2tau2 = 1;
    dim->v2rholapl = dim->v2rhotau = dim->v2lapltau = 1;
    dim->v2sigmalapl = dim->v2sigmatau = 1;

    dim->v3lapl3 = dim->v3tau3 = dim->v3rho2lapl = dim->v3rho2tau = dim->v3rholapl2 = 1;
    dim->v3rhotau2 = dim->v3lapl2tau = dim->v3lapltau2 = dim->v3rholapltau = 1;
    dim->v3sigmalapl2 = dim->v3sigmatau2 = dim->v3sigma2lapl = dim->v3sigma2tau = 1;
    dim->v3rhosigmalapl = dim->v3rhosigmatau = dim->v3sigmalapltau = 1;

    dim->v4rho4 = dim->v4rho3sigma = dim->v4rho3lapl = dim->v4rho3tau = dim->v4rho2sigma2 = 1;
    dim->v4rho2sigmalapl = dim->v4rho2sigmatau = dim->v4rho2lapl2 = dim->v4rho2lapltau = 1;
    dim->v4rho2tau2 = dim->v4rhosigma3 = dim->v4rhosigma2lapl = dim->v4rhosigma2tau = 1;
    dim->v4rhosigmalapl2 = dim->v4rhosigmalapltau = dim->v4rhosigmatau2 = 1;
    dim->v4rholapl3 = dim->v4rholapl2tau = dim->v4rholapltau2 = dim->v4rhotau3 = 1;
    dim->v4sigma4 = dim->v4sigma3lapl = dim->v4sigma3tau = dim->v4sigma2lapl2 = 1;
    dim->v4sigma2lapltau = dim->v4sigma2tau2 = dim->v4sigmalapl3 = dim->v4sigmalapl2tau = 1;
    dim->v4sigmalapltau2 = dim->v4sigmatau3 = dim->v4lapl4 = dim->v4lapl3tau = 1;
    dim->v4lapl2tau2 = dim->v4lapltau3 = dim->v4tau4 =1;
  }else{
    /* in total: 30 */
    dim->v2rholapl   = 2*2;
    dim->v2rhotau    = 2*2;
    dim->v2sigmalapl = 3*2;
    dim->v2sigmatau  = 3*2;
    dim->v2lapl2     = 3;
    dim->v2lapltau   = 2*2;
    dim->v2tau2      = 3;

    /* in total: 130 */
    dim->v3rho2lapl     = 3*2;
    dim->v3rho2tau      = 3*2;
    dim->v3rhosigmalapl = 2*3*2;
    dim->v3rhosigmatau  = 2*3*2;
    dim->v3rholapl2     = 2*3;
    dim->v3rholapltau   = 2*2*2;
    dim->v3rhotau2      = 2*3;
    dim->v3sigma2lapl   = 6*2;
    dim->v3sigma2tau    = 6*2;
    dim->v3sigmalapl2   = 3*3;
    dim->v3sigmalapltau = 3*2*2;
    dim->v3sigmatau2    = 3*3;
    dim->v3lapl3        = 4;
    dim->v3lapl2tau     = 3*2;
    dim->v3lapltau2     = 2*3;
    dim->v3tau3         = 4;

    /* in total: 477 */
    dim->v4rho3lapl        = 4*2;
    dim->v4rho3tau         = 4*2;
    dim->v4rho2sigmalapl   = 3*3*2;
    dim->v4rho2sigmatau    = 3*3*2;
    dim->v4rho2lapl2       = 3*3;
    dim->v4rho2lapltau     = 3*2*2;
    dim->v4rho2tau2        = 3*3;
    dim->v4rhosigma2lapl   = 3*6*2;
    dim->v4rhosigma2tau    = 3*6*2;
    dim->v4rhosigmalapl2   = 2*3*3;
    dim->v4rhosigmalapltau = 2*3*2*2;
    dim->v4rhosigmatau2    = 2*6*3;
    dim->v4rholapl3        = 2*4;
    dim->v4rholapl2tau     = 2*3*2;
    dim->v4rholapltau2     = 2*2*3;
    dim->v4rhotau3         = 2*4;
    dim->v4sigma3lapl      = 10*2;
    dim->v4sigma3tau       = 10*3;
    dim->v4sigma2lapl2     = 6*3;
    dim->v4sigma2lapltau   = 6*2*2;
    dim->v4sigma2tau2      = 6*3;
    dim->v4sigmalapl3      = 3*4;
    dim->v4sigmalapl2tau   = 3*3*2;
    dim->v4sigmalapltau2   = 3*2*3;
    dim->v4sigmatau3       = 3*4;
    dim->v4lapl4           = 5;
    dim->v4lapl3tau        = 4*2;
    dim->v4lapl2tau2       = 3*3;
    dim->v4lapltau3        = 2*4;
    dim->v4tau4            = 5;
  }
}

GPU_FUNCTION void
internal_counters_lda_random
  (const xc_dimensions *dim, int pos, int offset, const double **rho,
   double **zk, LDA_OUT_PARAMS_NO_EXC(double **))
{
  *rho += pos*dim->rho + offset;
  if(*zk != NULL)     *zk     += pos*dim->zk     + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL)   *vrho   += pos*dim->vrho   + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) *v2rho2 += pos*dim->v2rho2 + offset;
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) *v3rho3 += pos*dim->v3rho3 + offset;
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) *v4rho4 += pos*dim->v4rho4 + offset;
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_lda_next
  (const xc_dimensions *dim, int offset, const double **rho,
   double **zk, LDA_OUT_PARAMS_NO_EXC(double **))
{
  *rho += dim->rho + offset;
  if(*zk != NULL)     *zk     += dim->zk     + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL)   *vrho   += dim->vrho   + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) *v2rho2 += dim->v2rho2 + offset;
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) *v3rho3 += dim->v3rho3 + offset;
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) *v4rho4 += dim->v4rho4 + offset;
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_lda_prev
  (const xc_dimensions *dim, int offset, const double **rho,
   double **zk, LDA_OUT_PARAMS_NO_EXC(double **))
{
  *rho -= dim->rho + offset;
  if(*zk != NULL)     *zk     -= dim->zk     + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL)   *vrho   -= dim->vrho   + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) *v2rho2 -= dim->v2rho2 + offset;
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) *v3rho3 -= dim->v3rho3 + offset;
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) *v4rho4 -= dim->v4rho4 + offset;
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_gga_random
  (
   const xc_dimensions *dim, int pos, int offset, const double **rho, const double **sigma,
   double **zk, GGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_lda_random(dim, pos, offset, rho, zk, LDA_OUT_PARAMS_NO_EXC(XC_NOARG));

  *sigma += pos*dim->sigma + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) *vsigma += pos*dim->vsigma + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    *v2rhosigma += pos*dim->v2rhosigma + offset;
    *v2sigma2   += pos*dim->v2sigma2  + offset;
  }
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    *v3rho2sigma += pos*dim->v3rho2sigma + offset;
    *v3rhosigma2 += pos*dim->v3rhosigma2 + offset;
    *v3sigma3    += pos*dim->v3sigma3    + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    *v4rho3sigma  += pos*dim->v4rho3sigma  + offset;
    *v4rho2sigma2 += pos*dim->v4rho2sigma2 + offset;
    *v4rhosigma3  += pos*dim->v4rhosigma3  + offset;
    *v4sigma4     += pos*dim->v4sigma4     + offset;
  }
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_gga_next
  (
   const xc_dimensions *dim, int offset, const double **rho, const double **sigma,
   double **zk, GGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_lda_next(dim, offset, rho, zk, LDA_OUT_PARAMS_NO_EXC(XC_NOARG));

  *sigma += dim->sigma + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) *vsigma += dim->vsigma + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    *v2rhosigma += dim->v2rhosigma + offset;
    *v2sigma2   += dim->v2sigma2  + offset;
  }
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    *v3rho2sigma += dim->v3rho2sigma + offset;
    *v3rhosigma2 += dim->v3rhosigma2 + offset;
    *v3sigma3    += dim->v3sigma3    + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    *v4rho3sigma  += dim->v4rho3sigma  + offset;
    *v4rho2sigma2 += dim->v4rho2sigma2 + offset;
    *v4rhosigma3  += dim->v4rhosigma3  + offset;
    *v4sigma4     += dim->v4sigma4     + offset;
  }
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_gga_prev
(const xc_dimensions *dim, int offset, const double **rho, const double **sigma,
 double **zk, GGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_lda_prev(dim, offset, rho, zk, LDA_OUT_PARAMS_NO_EXC(XC_NOARG));

  *sigma -= dim->sigma + offset;
#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) *vsigma -= dim->vsigma   + offset;
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    *v2rhosigma -= dim->v2rhosigma + offset;
    *v2sigma2   -= dim->v2sigma2  + offset;
  }
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    *v3rho2sigma -= dim->v3rho2sigma + offset;
    *v3rhosigma2 -= dim->v3rhosigma2 + offset;
    *v3sigma3    -= dim->v3sigma3    + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    *v4rho3sigma  -= dim->v4rho3sigma  + offset;
    *v4rho2sigma2 -= dim->v4rho2sigma2 + offset;
    *v4rhosigma3  -= dim->v4rhosigma3  + offset;
    *v4sigma4     -= dim->v4sigma4     + offset;
  }
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_mgga_random
  (const xc_dimensions *dim, int pos, int offset,
   const double **rho, const double **sigma, const double **lapl, const double **tau,
   double **zk, MGGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_gga_random(dim, pos, offset, rho, sigma, zk, GGA_OUT_PARAMS_NO_EXC(XC_NOARG));

  if (*lapl != NULL)
    *lapl += pos*dim->lapl + offset;
  *tau  += pos*dim->tau  + offset;

#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) {
    if (*lapl != NULL)
      *vlapl += pos*dim->vlapl + offset;
    *vtau  += pos*dim->vtau  + offset;
  }
  
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    if (*lapl != NULL){
      *v2rholapl   += pos*dim->v2rholapl   + offset;
      *v2sigmalapl += pos*dim->v2sigmalapl + offset;
      *v2lapl2     += pos*dim->v2lapl2     + offset;
      *v2lapltau   += pos*dim->v2lapltau   + offset;
    }
    *v2rhotau    += pos*dim->v2rhotau    + offset;
    *v2sigmatau  += pos*dim->v2sigmatau  + offset;
    *v2tau2      += pos*dim->v2tau2      + offset;
  }
  
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    if (*lapl != NULL){
      *v3rho2lapl     += pos*dim->v3rho2lapl     + offset;
      *v3rhosigmalapl += pos*dim->v3rhosigmalapl + offset;
      *v3rholapl2     += pos*dim->v3rholapl2     + offset;
      *v3rholapltau   += pos*dim->v3rholapltau   + offset;
      *v3sigma2lapl   += pos*dim->v3sigma2lapl   + offset;
      *v3sigmalapl2   += pos*dim->v3sigmalapl2   + offset;
      *v3sigmalapltau += pos*dim->v3sigmalapltau + offset;
      *v3lapl3        += pos*dim->v3lapl3        + offset;
      *v3lapl2tau     += pos*dim->v3lapl2tau     + offset;
      *v3lapltau2     += pos*dim->v3lapltau2     + offset;
    }
    *v3rho2tau      += pos*dim->v3rho2tau      + offset;
    *v3rhosigmatau  += pos*dim->v3rhosigmatau  + offset;
    *v3rhotau2      += pos*dim->v3rhotau2      + offset;
    *v3sigma2tau    += pos*dim->v3sigma2tau    + offset;
    *v3sigmatau2    += pos*dim->v3sigmatau2    + offset;
    *v3tau3         += pos*dim->v3tau3         + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    if (*lapl != NULL){
      *v4rho3lapl        += pos*dim->v4rho3lapl        + offset;
      *v4rho2sigmalapl   += pos*dim->v4rho2sigmalapl   + offset;
      *v4rho2lapl2       += pos*dim->v4rho2lapl2       + offset;
      *v4rho2lapltau     += pos*dim->v4rho2lapltau     + offset;
      *v4rhosigma2lapl   += pos*dim->v4rhosigma2lapl   + offset;
      *v4rhosigmalapl2   += pos*dim->v4rhosigmalapl2   + offset;
      *v4rhosigmalapltau += pos*dim->v4rhosigmalapltau + offset;
      *v4rholapl3        += pos*dim->v4rholapl3        + offset;
      *v4rholapl2tau     += pos*dim->v4rholapl2tau     + offset;
      *v4rholapltau2     += pos*dim->v4rholapltau2     + offset;
      *v4sigma3lapl      += pos*dim->v4sigma3lapl      + offset;
      *v4sigma2lapl2     += pos*dim->v4sigma2lapl2     + offset;
      *v4sigma2lapltau   += pos*dim->v4sigma2lapltau   + offset;
      *v4sigmalapl3      += pos*dim->v4sigmalapl3      + offset;
      *v4sigmalapl2tau   += pos*dim->v4sigmalapl2tau   + offset;
      *v4sigmalapltau2   += pos*dim->v4sigmalapltau2   + offset;
      *v4lapl4           += pos*dim->v4lapl4           + offset;
      *v4lapl3tau        += pos*dim->v4lapl3tau        + offset;
      *v4lapl2tau2       += pos*dim->v4lapl2tau2       + offset;
      *v4lapltau3        += pos*dim->v4lapltau3        + offset;
    }
    *v4rho3tau         += pos*dim->v4rho3tau         + offset;
    *v4rho2sigmatau    += pos*dim->v4rho2sigmatau    + offset;
    *v4rho2tau2        += pos*dim->v4rho2tau2        + offset;
    *v4rhosigma2tau    += pos*dim->v4rhosigma2tau    + offset;
    *v4rhosigmatau2    += pos*dim->v4rhosigmatau2    + offset;
    *v4rhotau3         += pos*dim->v4rhotau3         + offset;
    *v4sigma3tau       += pos*dim->v4sigma3tau       + offset;
    *v4sigma2tau2      += pos*dim->v4sigma2tau2      + offset;
    *v4sigmatau3       += pos*dim->v4sigmatau3       + offset;
    *v4tau4            += pos*dim->v4tau4            + offset;
  }
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_mgga_next
  (const xc_dimensions *dim, int offset,
   const double **rho, const double **sigma, const double **lapl, const double **tau,
   double **zk, MGGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_gga_next(dim, offset, rho, sigma, zk, GGA_OUT_PARAMS_NO_EXC(XC_NOARG));

  if (*lapl != NULL)
    *lapl += dim->lapl + offset;
  *tau  += dim->tau  + offset;

#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) {
    if (*lapl != NULL)
      *vlapl += dim->vlapl + offset;
    *vtau  += dim->vtau  + offset;
  }
  
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    if (*lapl != NULL){
      *v2rholapl   += dim->v2rholapl   + offset;
      *v2sigmalapl += dim->v2sigmalapl + offset;
      *v2lapl2     += dim->v2lapl2     + offset;
      *v2lapltau   += dim->v2lapltau   + offset;
    }
    *v2rhotau    += dim->v2rhotau    + offset;
    *v2sigmatau  += dim->v2sigmatau  + offset;
    *v2tau2      += dim->v2tau2      + offset;
  }
  
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    if (*lapl != NULL){
      *v3rho2lapl     += dim->v3rho2lapl     + offset;
      *v3rhosigmalapl += dim->v3rhosigmalapl + offset;
      *v3rholapl2     += dim->v3rholapl2     + offset;
      *v3rholapltau   += dim->v3rholapltau   + offset;
      *v3sigma2lapl   += dim->v3sigma2lapl   + offset;
      *v3sigmalapl2   += dim->v3sigmalapl2   + offset;
      *v3sigmalapltau += dim->v3sigmalapltau + offset;
      *v3lapl3        += dim->v3lapl3        + offset;
      *v3lapl2tau     += dim->v3lapl2tau     + offset;
      *v3lapltau2     += dim->v3lapltau2     + offset;
    }
    *v3rho2tau      += dim->v3rho2tau      + offset;
    *v3rhosigmatau  += dim->v3rhosigmatau  + offset;
    *v3rhotau2      += dim->v3rhotau2      + offset;
    *v3sigma2tau    += dim->v3sigma2tau    + offset;
    *v3sigmatau2    += dim->v3sigmatau2    + offset;
    *v3tau3         += dim->v3tau3         + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    if (*lapl != NULL){
      *v4rho3lapl        += dim->v4rho3lapl        + offset;
      *v4rho2sigmalapl   += dim->v4rho2sigmalapl   + offset;
      *v4rho2lapl2       += dim->v4rho2lapl2       + offset;
      *v4rho2lapltau     += dim->v4rho2lapltau     + offset;
      *v4rhosigma2lapl   += dim->v4rhosigma2lapl   + offset;
      *v4rhosigmalapl2   += dim->v4rhosigmalapl2   + offset;
      *v4rhosigmalapltau += dim->v4rhosigmalapltau + offset;
      *v4rholapl3        += dim->v4rholapl3        + offset;
      *v4rholapl2tau     += dim->v4rholapl2tau     + offset;
      *v4rholapltau2     += dim->v4rholapltau2     + offset;
      *v4sigma3lapl      += dim->v4sigma3lapl      + offset;
      *v4sigma2lapl2     += dim->v4sigma2lapl2     + offset;
      *v4sigma2lapltau   += dim->v4sigma2lapltau   + offset;
      *v4sigmalapl3      += dim->v4sigmalapl3      + offset;
      *v4sigmalapl2tau   += dim->v4sigmalapl2tau   + offset;
      *v4sigmalapltau2   += dim->v4sigmalapltau2   + offset;
      *v4lapl4           += dim->v4lapl4           + offset;
      *v4lapl3tau        += dim->v4lapl3tau        + offset;
      *v4lapl2tau2       += dim->v4lapl2tau2       + offset;
      *v4lapltau3        += dim->v4lapltau3        + offset;
    }
    *v4rho3tau         += dim->v4rho3tau         + offset;
    *v4rho2sigmatau    += dim->v4rho2sigmatau    + offset;
    *v4rho2tau2        += dim->v4rho2tau2        + offset;
    *v4rhosigma2tau    += dim->v4rhosigma2tau    + offset;
    *v4rhosigmatau2    += dim->v4rhosigmatau2    + offset;
    *v4rhotau3         += dim->v4rhotau3         + offset;
    *v4sigma3tau       += dim->v4sigma3tau       + offset;
    *v4sigma2tau2      += dim->v4sigma2tau2      + offset;
    *v4sigmatau3       += dim->v4sigmatau3       + offset;
    *v4tau4            += dim->v4tau4            + offset;
  }
#endif
#endif
#endif
#endif
}

GPU_FUNCTION void
internal_counters_mgga_prev
  (const xc_dimensions *dim, int offset,
   const double **rho, const double **sigma, const double **lapl, const double **tau,
   double **zk, MGGA_OUT_PARAMS_NO_EXC(double **))
{
  internal_counters_gga_prev(dim, offset, rho, sigma, zk, GGA_OUT_PARAMS_NO_EXC(XC_NOARG));

  if(*lapl != NULL)
    *lapl -= dim->lapl + offset;
  *tau  -= dim->tau  + offset;

#ifndef XC_DONT_COMPILE_VXC
  if(*vrho != NULL) {
    if(*lapl != NULL)
      *vlapl -= dim->vlapl + offset;
    *vtau  -= dim->vtau  + offset;
  }
  
#ifndef XC_DONT_COMPILE_FXC
  if(*v2rho2 != NULL) {
    if(*lapl != NULL){
      *v2rholapl   -= dim->v2rholapl   + offset;
      *v2sigmalapl -= dim->v2sigmalapl + offset;
      *v2lapl2     -= dim->v2lapl2     + offset;
      *v2lapltau   -= dim->v2lapltau   + offset;
    }
    *v2rhotau    -= dim->v2rhotau    + offset;
    *v2sigmatau  -= dim->v2sigmatau  + offset;
    *v2tau2      -= dim->v2tau2      + offset;
  }
  
#ifndef XC_DONT_COMPILE_KXC
  if(*v3rho3 != NULL) {
    if (*lapl != NULL){
      *v3rho2lapl     -= dim->v3rho2lapl     + offset;
      *v3rhosigmalapl -= dim->v3rhosigmalapl + offset;
      *v3rholapl2     -= dim->v3rholapl2     + offset;
      *v3rholapltau   -= dim->v3rholapltau   + offset;
      *v3sigma2lapl   -= dim->v3sigma2lapl   + offset;
      *v3sigmalapl2   -= dim->v3sigmalapl2   + offset;
      *v3sigmalapltau -= dim->v3sigmalapltau + offset;
      *v3lapl3        -= dim->v3lapl3        + offset;
      *v3lapl2tau     -= dim->v3lapl2tau     + offset;
      *v3lapltau2     -= dim->v3lapltau2     + offset;
    }
    *v3rho2tau      -= dim->v3rho2tau      + offset;
    *v3rhosigmatau  -= dim->v3rhosigmatau  + offset;
    *v3rhotau2      -= dim->v3rhotau2      + offset;
    *v3sigma2tau    -= dim->v3sigma2tau    + offset;
    *v3sigmatau2    -= dim->v3sigmatau2    + offset;
    *v3tau3         -= dim->v3tau3         + offset;
  }
#ifndef XC_DONT_COMPILE_LXC
  if(*v4rho4 != NULL) {
    if (*lapl != NULL){
      *v4rho3lapl        -= dim->v4rho3lapl        + offset;
      *v4rho2sigmalapl   -= dim->v4rho2sigmalapl   + offset;
      *v4rho2lapl2       -= dim->v4rho2lapl2       + offset;
      *v4rho2lapltau     -= dim->v4rho2lapltau     + offset;
      *v4rhosigma2lapl   -= dim->v4rhosigma2lapl   + offset;
      *v4rhosigmalapl2   -= dim->v4rhosigmalapl2   + offset;
      *v4rhosigmalapltau -= dim->v4rhosigmalapltau + offset;
      *v4rholapl3        -= dim->v4rholapl3        + offset;
      *v4rholapl2tau     -= dim->v4rholapl2tau     + offset;
      *v4rholapltau2     -= dim->v4rholapltau2     + offset;
      *v4sigma3lapl      -= dim->v4sigma3lapl      + offset;
      *v4sigma2lapl2     -= dim->v4sigma2lapl2     + offset;
      *v4sigma2lapltau   -= dim->v4sigma2lapltau   + offset;
      *v4sigmalapl3      -= dim->v4sigmalapl3      + offset;
      *v4sigmalapl2tau   -= dim->v4sigmalapl2tau   + offset;
      *v4sigmalapltau2   -= dim->v4sigmalapltau2   + offset;
      *v4lapl4           -= dim->v4lapl4           + offset;
      *v4lapl3tau        -= dim->v4lapl3tau        + offset;
      *v4lapl2tau2       -= dim->v4lapl2tau2       + offset;
      *v4lapltau3        -= dim->v4lapltau3        + offset;
    }
    *v4rho3tau         -= dim->v4rho3tau         + offset;
    *v4rho2sigmatau    -= dim->v4rho2sigmatau    + offset;
    *v4rho2tau2        -= dim->v4rho2tau2        + offset;
    *v4rhosigma2tau    -= dim->v4rhosigma2tau    + offset;
    *v4rhosigmatau2    -= dim->v4rhosigmatau2    + offset;
    *v4rhotau3         -= dim->v4rhotau3         + offset;
    *v4sigma3tau       -= dim->v4sigma3tau       + offset;
    *v4sigma2tau2      -= dim->v4sigma2tau2      + offset;
    *v4sigmatau3       -= dim->v4sigmatau3       + offset;
    *v4tau4            -= dim->v4tau4            + offset;
  }
#endif
#endif
#endif
#endif
}
