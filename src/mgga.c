/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_mgga.c"

#define FUNC(type)   xc_mgga ## type

#define IN_VARIABLES                         \
  double *rho, double *sigma, double *lapl, double *tau
#define INPUT_VARIABLES                      \
  rho, sigma, lapl, tau, NULL


#define OUT_VARIABLES_0                      \
  double *zk
#define SET_ORDER_0                          \
  out.zk = zk;


#define OUT_VARIABLES_1                      \
  double *vrho, double *vsigma,              \
  double *vlapl, double *vtau
#define SET_ORDER_1                          \
  out.vrho = vrho;                           \
  out.vsigma = vsigma;                       \
  out.vlapl = vlapl;                         \
  out.vtau = vtau;


#define OUT_VARIABLES_2                      \
  double *v2rho2, double *v2rhosigma,        \
  double *v2rholapl, double *v2rhotau,       \
  double *v2sigma2, double *v2sigmalapl,     \
  double *v2sigmatau, double *v2lapl2,       \
  double *v2lapltau, double *v2tau2
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


#define OUT_VARIABLES_3                      \
   double *v3rho3, double *v3rho2sigma,      \
   double *v3rho2lapl, double *v3rho2tau,    \
   double *v3rhosigma2, double *v3rhosigmalapl, \
   double *v3rhosigmatau, double *v3rholapl2,   \
   double *v3rholapltau, double *v3rhotau2,     \
   double *v3sigma3, double *v3sigma2lapl,      \
   double *v3sigma2tau, double *v3sigmalapl2,   \
   double *v3sigmalapltau, double *v3sigmatau2, \
   double *v3lapl3, double *v3lapl2tau,         \
   double *v3lapltau2, double *v3tau3
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


#define OUT_VARIABLES_4                              \
  double *v4rho4, double *v4rho3sigma,               \
  double *v4rho3lapl, double *v4rho3tau,             \
  double *v4rho2sigma2, double *v4rho2sigmalapl,     \
  double *v4rho2sigmatau, double *v4rho2lapl2,       \
  double *v4rho2lapltau, double *v4rho2tau2,         \
  double *v4rhosigma3, double *v4rhosigma2lapl,      \
  double *v4rhosigma2tau, double *v4rhosigmalapl2,   \
  double *v4rhosigmalapltau, double *v4rhosigmatau2, \
  double *v4rholapl3, double *v4rholapl2tau,         \
  double *v4rholapltau2, double *v4rhotau3,          \
  double *v4sigma4, double *v4sigma3lapl,            \
  double *v4sigma3tau, double *v4sigma2lapl2,        \
  double *v4sigma2lapltau, double *v4sigma2tau2,     \
  double *v4sigmalapl3, double *v4sigmalapl2tau,     \
  double *v4sigmalapltau2, double *v4sigmatau3,      \
  double *v4lapl4, double *v4lapl3tau,               \
  double *v4lapl2tau2, double *v4lapltau3,           \
  double *v4tau4
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


#include "old_interface.c"
