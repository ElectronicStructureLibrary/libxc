/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_gga.c"

#define FUNC(type)   xc_gga ## type

#define IN_VARIABLES                         \
  double *rho, double *sigma
#define INPUT_VARIABLES                      \
  rho, sigma, NULL, NULL, NULL


#define OUT_VARIABLES_0                      \
  double *zk
#define SET_ORDER_0                          \
  out.zk = zk;


#define OUT_VARIABLES_1                      \
  double *vrho, double *vsigma
#define SET_ORDER_1                          \
  out.vrho = vrho;                           \
  out.vsigma = vsigma;

#define OUT_VARIABLES_2                      \
  double *v2rho2, double *v2rhosigma, double *v2sigma2
#define SET_ORDER_2                          \
  out.v2rho2 = v2rho2;                       \
  out.v2rhosigma = v2rhosigma;               \
  out.v2sigma2 = v2sigma2;

#define OUT_VARIABLES_3                      \
  double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3
#define SET_ORDER_3                          \
  out.v3rho3 = v3rho3;                       \
  out.v3rho2sigma = v3rho2sigma;             \
  out.v3rhosigma2 = v3rhosigma2;             \
  out.v3sigma3 = v3sigma3;

#define OUT_VARIABLES_4                      \
  double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4
#define SET_ORDER_4                          \
  out.v4rho4 = v4rho4;                       \
  out.v4rho3sigma = v4rho3sigma;             \
  out.v4rho2sigma2 = v4rho2sigma2;           \
  out.v4rhosigma3 = v4rhosigma3;             \
  out.v4sigma4 = v4sigma4;


#include "old_interface.c"
