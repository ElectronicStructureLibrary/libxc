/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include "util.h"

/************************************************************************
 Random Phase Approximation (RPA)
************************************************************************/

#define XC_LDA_C_RPA  3   /* Random Phase Approximation   */

static inline void 
func(const XC(lda_type) *p, int order, FLOAT *rs, FLOAT zeta, 
     FLOAT *zk, FLOAT *dedrs, FLOAT *dedz, 
     FLOAT *d2edrs2, FLOAT *d2edrsz, FLOAT *d2edz2)
{
  static FLOAT a = 0.0311, b = -0.047, c = 0.009, d = -0.017;
  FLOAT lrs;

  lrs = log(rs[1]);
  *zk = a*lrs + b + c*rs[1]*lrs + d*rs[1];

  if(order < 1) return;

  *dedrs = a/rs[1] + c*(lrs + 1.0) + d;

  if(order < 2) return;

  *d2edrs2 = -a/rs[2] + c/rs[1];
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_rpa) = {
  XC_LDA_C_RPA,
  XC_CORRELATION,
  "Random Phase Approximation (RPA)",
  XC_FAMILY_LDA,
  "M Gell-Mann and KA Brueckner, Phys. Rev. 106, 364 (1957)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
