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

#define XC_HYB_GGA_XC_PBEH 406 /* aka PBE0 or PBE1PBE */

static void
gga_xc_pbeh_init(void *p_)
{
  const FLOAT a0 = 0.25;

  XC(hyb_gga_type) *p = (XC(hyb_gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 0;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  p->exx_coef = a0;

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_PBE, p->nspin);
  p->mix->gga_coef[0] = (1.0 - a0);
  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_C_PBE, p->nspin);
  p->mix->gga_coef[1] = 1.0;
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_pbeh) = {
  XC_HYB_GGA_XC_PBEH,
  XC_EXCHANGE_CORRELATION,
  "PBEH (PBE0)",
  XC_FAMILY_HYB_GGA,
  "M. Ernzerhof, G. E. Scuseria, J. Chem. Phys. 110, 5029 (1999)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_pbeh_init,
  NULL, 
  NULL,
  NULL /* this is taken care by the generic routine */
};
