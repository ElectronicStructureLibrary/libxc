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

#define XC_GGA_XC_XLYP 166 /* XLYP functional */

static void
gga_xc_xlyp_init(void *p_)
{
  const FLOAT cc[3] = {0.722, 0.347, 1.0};

  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 1;
  p->mix->gga_n = 3;
  XC(mix_func_alloc)(p->mix);

  XC(lda_init)(&p->mix->lda_mix[0], XC_LDA_X, p->nspin);
  p->mix->lda_coef[0] = 1.0 - cc[0] - cc[1];

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_B88,  p->nspin);
  p->mix->gga_coef[0] = cc[0];

  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_X_PW91, p->nspin);
  p->mix->gga_coef[0] = cc[1];

  XC(gga_init)(&p->mix->gga_mix[2], XC_GGA_C_LYP,  p->nspin);
  p->mix->gga_coef[0] = cc[2];
}

const XC(func_info_type) XC(func_info_gga_xc_xlyp) = {
  XC_GGA_XC_XLYP,
  XC_EXCHANGE_CORRELATION,
  "XLYP",
  XC_FAMILY_GGA,
  "X Xu and WA Goddard, III, PNAS 101, 2673 (2004)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_xlyp_init, 
  NULL, NULL, NULL
};
