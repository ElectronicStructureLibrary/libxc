/*
 Copyright (C) 2006-2016 M.A.L. Marques

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


#include "util.h"

#define XC_LDA_C_CHACHIYO  287   /* Chachiyo simple 2 parameter correlation   */

void
XC(lda_c_chachiyo_func)(const XC(func_type) *p, XC(lda_work_t) *r)
{
  static FLOAT ap = -0.01554535, bp = 20.4562557;
  static FLOAT af = -.007772675, bf = 27.4203609; /* af = ap/2.0 */

  FLOAT aux0, daux0, d2aux0, d3aux0, e0, de0, d2e0, d3e0;
  FLOAT aux1, daux1, d2aux1, d3aux1, e1, de1, d2e1, d3e1;
  FLOAT f, dfdz, d2fdz2, d3fdz3;

  /* paramagnetic case */
  aux0 = 1.0 + bp/r->rs[1] + bp/r->rs[2];
  e0 = ap*LOG(aux0);

  /* ferromagnetic case (hartree) */
  aux1 = 1.0 + bf/r->rs[1] + bf/r->rs[2];
  e1 = af*LOG(aux1);

  /* spin polarization */
  f = FZETA(r->zeta);

  r->zk  = e0 + (e1 - e0)*f;

  if(r->order < 1) return;

  daux0 = -bp/r->rs[2] - 2.0*bp/(r->rs[1]*r->rs[2]);
  de0 = ap*daux0/aux0;

  daux1 = -bf/r->rs[2] - 2.0*bf/(r->rs[1]*r->rs[2]);
  de1 = af*daux1/aux1;
  
  dfdz = DFZETA(r->zeta);
  
  r->dedrs = de0 + (de1 - de0)*f;
  r->dedz = (e1 - e0)*dfdz;

  if(r->order < 2) return;

  d2aux0 = 2.0*bp/(r->rs[1]*r->rs[2]) + 6.0*bp/(r->rs[2]*r->rs[2]);
  d2e0   = ap*DFRACTION(daux0, d2aux0, aux0, daux0);

  d2aux1 = 2.0*bf/(r->rs[1]*r->rs[2]) + 6.0*bf/(r->rs[2]*r->rs[2]);
  d2e1   = af*DFRACTION(daux1, d2aux1, aux1, daux1);

  d2fdz2 = D2FZETA(r->zeta);

  r->d2edrs2 = d2e0 + (d2e1 - d2e0)*f;
  r->d2edrsz = (de1 - de0)*dfdz;
  r->d2edz2  = (e1 - e0)*d2fdz2;

  if(r->order < 3) return;

  d3aux0 = -6.0*bp/(r->rs[2]*r->rs[2]) - 24.0*bp/(r->rs[1]*r->rs[2]*r->rs[2]);
  d3e0   = ap*D2FRACTION(daux0, d2aux0, d3aux0, aux0, daux0, d2aux0);

  d3aux1 = -6.0*bf/(r->rs[2]*r->rs[2]) - 24.0*bf/(r->rs[1]*r->rs[2]*r->rs[2]);
  d3e1   = af*D2FRACTION(daux1, d2aux1, d3aux1, aux1, daux1, d2aux1);

  d3fdz3 = D3FZETA(r->zeta);

  r->d3edrs3  = d3e0 + (d3e1 - d3e0)*f;
  r->d3edrs2z = (d2e1 - d2e0)*dfdz;
  r->d3edrsz2 = (de1 - de0)*d2fdz2;
  r->d3edz3   = (e1 - e0)*d3fdz3;
}

#define func XC(lda_c_chachiyo_func)
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_chachiyo) = {
  XC_LDA_C_CHACHIYO,
  XC_CORRELATION,
  "Chachiyo simple 2 parameter correlation",
  XC_FAMILY_LDA,
  {&xc_ref_Chachiyo2016_021101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
  NULL,
  NULL
};
