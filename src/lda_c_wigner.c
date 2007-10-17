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

static void lda_c_wigner(const void *p_, const double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a = -0.44, b = 7.8;
  double dens, zeta, rs;
  double etmp, decdrs, t;
  
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs    =  RS(dens); /* Wigner radius */
  t     =  b + rs;

  etmp   =  a/t;
  decdrs = -a/(t*t);                         /* now contains d ec/d rs */
  
  if(ec != NULL) *ec = etmp;

  if(vc != NULL){
    vc[0] = etmp - decdrs*rs/3.0;              /* and now d ec/d rho */
    if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to return something */
  }

}

const xc_func_info_type func_info_lda_c_wigner = {
  XC_LDA_C_WIGNER,
  XC_CORRELATION,
  "Wigner",
  XC_FAMILY_LDA,
  "E.P. Wigner, Trans. Faraday Soc. 34, 678 (1938)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_wigner, /* lda  */
};
