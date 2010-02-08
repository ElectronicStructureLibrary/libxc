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

#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Correlation energy of Proynov and Salahub
************************************************************************/

#define XC_LDA_C_PS94    22   /* Perdew & Zunger              */

/* the functional */
static inline void 
func(const XC(lda_type) *p, XC(lda_rs_zeta) *r)
{
  static FLOAT fc = 0.2026, q = 0.084, C = 6.187335;
  static FLOAT b[6] = {2.763169, 1.757515, 1.741397, 0.568985, 1.572202, 1.885389};

  FLOAT cnst_rs, zp3, zm3, alpha, k, Q;
  FLOAT dQ, dkdrs

  cnst_rs = POW(3.0/(4*M_PI), 1.0/3.0);

  alpha = fc*(pow(1 + r->zeta, q) + pow(1.0 - r->zeta, q));

  zp3   = pow(1.0 + r->zeta,  1.0/3.0);
  zm3   = pow(1.0 - r->zeta,  1.0/3.0);
  k     = C*alpha*cnst_rs/r->rs[1]*(zp3*zm3/(zp3 + zm3));

  Q = (k == 0.0) ? 0.0 : -b[0]/(1.0 + b[1]*k) + b[2]/k*log(1.0 + b[3]/k) + b[4]/k - b[5]/(k*k);

  r->zk = pow(cnst_rs/r->rs[1], 3) * (1 - r->zeta*r->zeta)/4.0 * Q;

  if(r->order < 1) return;

  dQ = b[0]*b[1]/((1.0 + b[1]*k)*(1.0 + b[1]*k)) + b[2]*b[3]/((b[3] + k)*(k*k))
    - b[2]*log(1.0 + b[3]/k)/(k*k) - b[4]/(k*k) + 2.0*b[5]/(k*k*k);
  dkdrs = -k/r->rs[1];

  r->dedrs = pow(cnst_rs/r->rs[1], 3)*(1 - r->zeta*r->zeta)/4.0 * (dQ*dkdrs + 3.0*Q/r->rs[1]);
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_ps94) = {
  XC_LDA_C_PS94,
  XC_CORRELATION,
  "Proynov and Salahub 94",
  XC_FAMILY_LDA,
  "EI Proynov and D Salahub, Phys. Rev. B 49, 7874 (1994)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
