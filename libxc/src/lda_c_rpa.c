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

static void lda_c_rpa(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static FLOAT a = 0.0311, b = -0.047, c = 0.009, d = -0.017;
  FLOAT dens, zeta, rs;
  FLOAT lrs;

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs  =  RS(dens); /* Wigner radius */
  lrs = log(rs);
  
  *ec   = a*lrs + b + c*rs*lrs + d*rs;
  vc[0] = a/rs + c*(lrs + 1.0) + d;         /* now contains d ec/d rs */
  
  vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
  if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}

const xc_func_info_type func_info_lda_c_rpa = {
  XC_LDA_C_RPA,
  XC_CORRELATION,
  "Random Phase Approximation (RPA)",
  XC_FAMILY_LDA,
  "M Gell-Mann and KA Brueckner, Phys. Rev. 106, 364 (1957)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_rpa,    /* lda  */
};
