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
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_HYB_GGA_XC_B3LYP 402 /* The (in)famous B3LYP */

static void
gga_xc_b3lyp_init(void *p_)
{
  const FLOAT a0 = 0.20, ax = 0.72, ac = 0.81;

  XC(hyb_gga_type) *p = (XC(hyb_gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 2;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  p->exx_coef = a0;

  XC(lda_x_init)(&p->mix->lda_mix[0], p->nspin, 3, XC_NON_RELATIVISTIC);
  p->mix->lda_coef[0] = 1.0 - a0 - ax;

  /* set the vwn part with the spin interpolation scheme originally used in Gaussian */
  XC(lda_init)  (&p->mix->lda_mix[1], XC_LDA_C_VWN_RPA, p->nspin);
  XC(lda_c_vwn_set_params)(&p->mix->lda_mix[1], 1);
  p->mix->lda_coef[1] = 1.0 - ac;

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_B88, p->nspin);
  p->mix->gga_coef[0] = ax;
  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_C_LYP, p->nspin);
  p->mix->gga_coef[1] = ac;
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyp) = {
  XC_HYB_GGA_XC_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "B3LYP",
  XC_FAMILY_HYB_GGA,
  "PJ Stephens, FJ Devlin, CF Chabalowski, MJ Frisch, J. Phys. Chem. 98 11623 (1994)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3lyp_init,
  NULL, NULL, NULL
};
