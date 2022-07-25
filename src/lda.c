/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_lda.c"

#define FUNC(type)   xc_lda ## type

#define IN_VARIABLES                         \
  double *rho
#define INPUT_VARIABLES                      \
  rho, NULL, NULL, NULL, NULL


#define OUT_VARIABLES_0                      \
  double *zk
#define SET_ORDER_0                          \
  out.zk = zk;


#define OUT_VARIABLES_1                      \
  double *vrho
#define SET_ORDER_1                          \
  out.vrho = vrho;


#define OUT_VARIABLES_2                      \
  double *v2rho2
#define SET_ORDER_2                          \
  out.v2rho2 = v2rho2;

#define OUT_VARIABLES_3                      \
  double *v3rho3
#define SET_ORDER_3                          \
  out.v3rho3 = v3rho3;

#define OUT_VARIABLES_4                      \
  double *v4rho4
#define SET_ORDER_4                          \
  out.v4rho4 = v4rho4;


#include "old_interface.c"
