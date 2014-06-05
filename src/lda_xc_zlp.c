/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_LDA_XC_ZLP     43   /* Zhao, Levy & Parr  */

/* the functional */
static inline void 
func(const XC(func_type) *p, XC(lda_work_t) *r)
{
  static FLOAT a0 = 0.93222, kk = 9.47362e-3/RS_FACTOR;
  FLOAT aux, daux, d2aux, d3aux;

  aux = LOG(1.0 + r->rs[1]/kk);

  r->zk = -a0*(1.0 - kk*aux/r->rs[1])*RS_FACTOR/r->rs[1];

  if(r->order < 1) return; /* nothing else to do */

  daux = 1.0/(r->rs[1] + kk);

  r->dedrs = a0*(r->rs[1] - 2.0*kk*aux + kk*r->rs[1]*daux)*RS_FACTOR/(r->rs[2]*r->rs[1]);
  r->dedz  = 0.0;

  if(r->order < 2) return; /* nothing else to do */

  d2aux = -daux*daux;

  r->d2edrs2  = a0*(-2.0*r->rs[1] + 6.0*kk*aux - 4.0*kk*r->rs[1]*daux + kk*r->rs[2]*d2aux)*RS_FACTOR/(r->rs[2]*r->rs[2]);
  r->d2edz2   = 0.0;
  r->d2edrsz  = 0.0;

  if(r->order < 3) return; /* nothing else to do */

  d3aux = -2.0*d2aux*daux;

  r->d3edrs3   = a0*(6.0*r->rs[1] - 24.0*kk*aux + 18.0*kk*r->rs[1]*daux
		     - 6.0*kk*r->rs[2]*d2aux + kk*r->rs[2]*r->rs[1]*d3aux)*RS_FACTOR/(r->rs[2]*r->rs[2]*r->rs[1]);
  r->d3edz3    = 0.0;
  r->d3edrs2z  = 0.0;
  r->d3edrsz2  = 0.0;
  r->d3edrsz2  = 0.0;
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_xc_zlp) = {
  XC_LDA_XC_ZLP,
  XC_EXCHANGE_CORRELATION,
  "Zhao, Levy & Parr",
  XC_FAMILY_LDA,
  {&xc_ref_Zhao1993_918, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
  NULL,
  NULL
};
