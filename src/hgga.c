/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_hgga.c"

#define FUNC(type)   xc_hgga ## type

#define IN_VARIABLES                         \
  double *rho, double *sigma, double *lapl, double *tau, double *exx
#define INPUT_VARIABLES                      \
  rho, sigma, lapl, tau, exx


#define OUT_VARIABLES_0                      \
  double *zk
#define SET_ORDER_0                          \
  out.zk = zk;


#define OUT_VARIABLES_1                      \
  double *vrho, double *vsigma,              \
  double *vlapl, double *vtau,               \
  double *vexx
#define SET_ORDER_1                          \
  out.vrho = vrho;                           \
  out.vsigma = vsigma;                       \
  out.vlapl = vlapl;                         \
  out.vtau = vtau;                           \
  out.vexx = vexx;


#define OUT_VARIABLES_2                      \
  double *v2rho2, double *v2rhosigma,        \
  double *v2rholapl, double *v2rhotau,       \
  double *v2rhoexx, double *v2sigma2,        \
  double *v2sigmalapl, double *v2sigmatau,   \
  double *v2sigmaexx, double *v2lapl2,       \
  double *v2lapltau, double *v2laplexx,      \
  double *v2tau2, double *v2tauexx,          \
  double *v2exx2
#define SET_ORDER_2                          \
  out.v2rho2 = v2rho2;                       \
  out.v2rhosigma = v2rhosigma;               \
  out.v2rholapl = v2rholapl;                 \
  out.v2rhotau = v2rhotau;                   \
  out.v2rhoexx = v2rhoexx;                   \
  out.v2sigma2 = v2sigma2;                   \
  out.v2sigmalapl = v2sigmalapl;             \
  out.v2sigmatau = v2sigmatau;               \
  out.v2sigmaexx = v2sigmaexx;               \
  out.v2lapl2 = v2lapl2;                     \
  out.v2lapltau = v2lapltau;                 \
  out.v2laplexx = v2laplexx;                 \
  out.v2tau2 = v2tau2;                       \
  out.v2tauexx = v2tauexx;                   \
  out.v2exx2 = v2exx2;


#define OUT_VARIABLES_3                             \
  double *v3rho3, double *v3rho2sigma,              \
  double *v3rho2lapl, double *v3rho2tau,            \
  double *v3rho2exx, double *v3rhosigma2,           \
  double *v3rhosigmalapl, double *v3rhosigmatau,    \
  double *v3rhosigmaexx, double *v3rholapl2,        \
  double *v3rholapltau, double *v3rholaplexx,       \
  double *v3rhotau2, double *v3rhotauexx,           \
  double *v3rhoexx2, double *v3sigma3,              \
  double *v3sigma2lapl, double *v3sigma2tau,        \
  double *v3sigma2exx, double *v3sigmalapl2,        \
  double *v3sigmalapltau, double *v3sigmalaplexx,   \
  double *v3sigmatau2, double *v3sigmatauexx,       \
  double *v3sigmaexx2, double *v3lapl3,             \
  double *v3lapl2tau, double *v3lapl2exx,           \
  double *v3lapltau2, double *v3lapltauexx,         \
  double *v3laplexx2, double *v3tau3,               \
  double *v3tau2exx, double *v3tauexx2,             \
  double *v3exx3
#define SET_ORDER_3                                 \
  out.v3rho3 = v3rho3;                              \
  out.v3rho2sigma = v3rho2sigma;                    \
  out.v3rho2lapl = v3rho2lapl;                      \
  out.v3rho2tau = v3rho2tau;                        \
  out.v3rho2exx = v3rho2exx;                        \
  out.v3rhosigma2 = v3rhosigma2;                    \
  out.v3rhosigmalapl = v3rhosigmalapl;              \
  out.v3rhosigmatau = v3rhosigmatau;                \
  out.v3rhosigmaexx = v3rhosigmaexx;                \
  out.v3rholapl2 = v3rholapl2;                      \
  out.v3rholapltau = v3rholapltau;                  \
  out.v3rholaplexx = v3rholaplexx;                  \
  out.v3rhotau2 = v3rhotau2;                        \
  out.v3rhotauexx = v3rhotauexx;                    \
  out.v3rhoexx2 = v3rhoexx2;                        \
  out.v3sigma3 = v3sigma3;                          \
  out.v3sigma2lapl = v3sigma2lapl;                  \
  out.v3sigma2tau = v3sigma2tau;                    \
  out.v3sigma2exx = v3sigma2exx;                    \
  out.v3sigmalapl2 = v3sigmalapl2;                  \
  out.v3sigmalapltau = v3sigmalapltau;              \
  out.v3sigmalaplexx = v3sigmalaplexx;              \
  out.v3sigmatau2 = v3sigmatau2;                    \
  out.v3sigmatauexx = v3sigmatauexx;                \
  out.v3sigmaexx2 = v3sigmaexx2;                    \
  out.v3lapl3 = v3lapl3;                            \
  out.v3lapl2tau = v3lapl2tau;                      \
  out.v3lapl2exx = v3lapl2exx;                      \
  out.v3lapltau2 = v3lapltau2;                      \
  out.v3lapltauexx = v3lapltauexx;                  \
  out.v3laplexx2 = v3laplexx2;                      \
  out.v3tau3 = v3tau3;                              \
  out.v3tau2exx = v3tau2exx;                        \
  out.v3tauexx2 = v3tauexx2;                        \
  out.v3exx3 = v3exx3;


#define OUT_VARIABLES_4                                   \
    double *v4rho4, double *v4rho3sigma,                  \
    double *v4rho3lapl, double *v4rho3tau,                \
    double *v4rho3exx, double *v4rho2sigma2,              \
    double *v4rho2sigmalapl, double *v4rho2sigmatau,      \
    double *v4rho2sigmaexx, double *v4rho2lapl2,          \
    double *v4rho2lapltau, double *v4rho2laplexx,         \
    double *v4rho2tau2, double *v4rho2tauexx,             \
    double *v4rho2exx2, double *v4rhosigma3,              \
    double *v4rhosigma2lapl, double *v4rhosigma2tau,      \
    double *v4rhosigma2exx, double *v4rhosigmalapl2,      \
    double *v4rhosigmalapltau, double *v4rhosigmalaplexx, \
    double *v4rhosigmatau2, double *v4rhosigmatauexx,     \
    double *v4rhosigmaexx2, double *v4rholapl3,           \
    double *v4rholapl2tau, double *v4rholapl2exx,         \
    double *v4rholapltau2, double *v4rholapltauexx,       \
    double *v4rholaplexx2, double *v4rhotau3,             \
    double *v4rhotau2exx, double *v4rhoexx3,              \
    double *v4sigma4, double *v4sigma3lapl,               \
    double *v4sigma3tau, double *v4sigma3exx,             \
    double *v4sigma2lapl2, double *v4sigma2lapltau,       \
    double *v4sigma2laplexx, double *v4sigma2tau2,        \
    double *v4sigma2tauexx, double *v4sigma2exx2,         \
    double *v4sigmalapl3, double *v4sigmalapl2tau,        \
    double *v4sigmalapl2exx, double *v4sigmalapltau2,     \
    double *v4sigmalapltauexx, double *v4sigmalaplexx2,   \
    double *v4sigmatau3, double *v4sigmatau2exx,          \
    double *v4sigmatauexx2, double *v4sigmaexx3,          \
    double *v4lapl4, double *v4lapl3tau,                  \
    double *v4lapl3exx, double *v4lapl2tau2,              \
    double *v4lapl2tauexx, double *v4lapl2exx2,           \
    double *v4lapltau3, double *v4lapltau2exx,            \
    double *v4lapltauexx2, double *v4laplexx3,            \
    double *v4tau4, double *v4tau3exx,                    \
    double *v4tauexx3, double *v4exx4
#define SET_ORDER_4                                       \
  out.v4rho4 = v4rho4;                                    \
  out.v4rho3sigma = v4rho3sigma;                          \
  out.v4rho3lapl = v4rho3lapl;                            \
  out.v4rho3tau = v4rho3tau;                              \
  out.v4rho3exx = v4rho3exx;                              \
  out.v4rho2sigma2 = v4rho2sigma2;                        \
  out.v4rho2sigmalapl = v4rho2sigmalapl;                  \
  out.v4rho2sigmatau = v4rho2sigmatau;                    \
  out.v4rho2sigmaexx = v4rho2sigmaexx;                    \
  out.v4rho2lapl2 = v4rho2lapl2;                          \
  out.v4rho2lapltau = v4rho2lapltau;                      \
  out.v4rho2laplexx = v4rho2laplexx;                      \
  out.v4rho2tau2 = v4rho2tau2;                            \
  out.v4rho2tauexx = v4rho2tauexx;                        \
  out.v4rho2exx2 = v4rho2exx2;                            \
  out.v4rhosigma3 = v4rhosigma3;                          \
  out.v4rhosigma2lapl = v4rhosigma2lapl;                  \
  out.v4rhosigma2tau = v4rhosigma2tau;                    \
  out.v4rhosigma2exx = v4rhosigma2exx;                    \
  out.v4rhosigmalapl2 = v4rhosigmalapl2;                  \
  out.v4rhosigmalapltau = v4rhosigmalapltau;              \
  out.v4rhosigmalaplexx = v4rhosigmalaplexx;              \
  out.v4rhosigmatau2 = v4rhosigmatau2;                    \
  out.v4rhosigmatauexx = v4rhosigmatauexx;                \
  out.v4rhosigmaexx2 = v4rhosigmaexx2;                    \
  out.v4rholapl3 = v4rholapl3;                            \
  out.v4rholapl2tau = v4rholapl2tau;                      \
  out.v4rholapl2exx = v4rholapl2exx;                      \
  out.v4rholapltau2 = v4rholapltau2;                      \
  out.v4rholapltauexx = v4rholapltauexx;                  \
  out.v4rholaplexx2 = v4rholaplexx2;                      \
  out.v4rhotau3 = v4rhotau3;                              \
  out.v4rhotau2exx = v4rhotau2exx;                        \
  out.v4rhoexx3 = v4rhoexx3;                              \
  out.v4sigma4 = v4sigma4;                                \
  out.v4sigma3lapl = v4sigma3lapl;                        \
  out.v4sigma3tau = v4sigma3tau;                          \
  out.v4sigma3exx = v4sigma3exx;                          \
  out.v4sigma2lapl2 = v4sigma2lapl2;                      \
  out.v4sigma2lapltau = v4sigma2lapltau;                  \
  out.v4sigma2laplexx = v4sigma2laplexx;                  \
  out.v4sigma2tau2 = v4sigma2tau2;                        \
  out.v4sigma2tauexx = v4sigma2tauexx;                    \
  out.v4sigma2exx2 = v4sigma2exx2;                        \
  out.v4sigmalapl3 = v4sigmalapl3;                        \
  out.v4sigmalapl2tau = v4sigmalapl2tau;                  \
  out.v4sigmalapl2exx = v4sigmalapl2exx;                  \
  out.v4sigmalapltau2 = v4sigmalapltau2;                  \
  out.v4sigmalapltauexx = v4sigmalapltauexx;              \
  out.v4sigmalaplexx2 = v4sigmalaplexx2;                  \
  out.v4sigmatau3 = v4sigmatau3;                          \
  out.v4sigmatau2exx = v4sigmatau2exx;                    \
  out.v4sigmatauexx2 = v4sigmatauexx2;                    \
  out.v4sigmaexx3 = v4sigmaexx3;                          \
  out.v4lapl4 = v4lapl4;                                  \
  out.v4lapl3tau = v4lapl3tau;                            \
  out.v4lapl3exx = v4lapl3exx;                            \
  out.v4lapl2tau2 = v4lapl2tau2;                          \
  out.v4lapl2tauexx = v4lapl2tauexx;                      \
  out.v4lapl2exx2 = v4lapl2exx2;                          \
  out.v4lapltau3 = v4lapltau3;                            \
  out.v4lapltau2exx = v4lapltau2exx;                      \
  out.v4lapltauexx2 = v4lapltauexx2;                      \
  out.v4laplexx3 = v4laplexx3;                            \
  out.v4tau4 = v4tau4;                                    \
  out.v4tau3exx = v4tau3exx;                              \
  out.v4tauexx3 = v4tauexx3;                              \
  out.v4exx4 = v4exx4;

#include "old_interface.c"
