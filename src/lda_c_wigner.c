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
 Wigner's parametrization from the low density limit
************************************************************************/

#define XC_LDA_C_WIGNER  2   /* Wigner parametrization       */

static inline void 
func(const XC(lda_type) *p, FLOAT *rs, FLOAT zeta, 
     FLOAT *zk, FLOAT *dedrs, FLOAT *dedz, 
     FLOAT *d2edrs2, FLOAT *d2edrsz, FLOAT *d2edz2)
{
  static FLOAT a = -0.44, b = 7.8;
  FLOAT t;
  
  t   =  b + rs[1];
  *zk =  a/t;

  if(dedrs != NULL)
    *dedrs = -a/(t*t);

  if(d2edrs2 != NULL)
    *d2edrs2 = 2.0*a/(t*t*t);
  
  if(p->nspin==XC_POLARIZED){
    if(dedrs   != NULL) *dedz = 0.0;
    if(d2edrs2 != NULL) *d2edrsz = *d2edz2 = 0.0;
  }

}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_wigner) = {
  XC_LDA_C_WIGNER,
  XC_CORRELATION,
  "Wigner",
  XC_FAMILY_LDA,
  "EP Wigner, Trans. Faraday Soc. 34, 678 (1938)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
