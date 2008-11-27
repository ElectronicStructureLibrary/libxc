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

#define XC_GGA_XC_EDF1 165 /* Empirical functionals from Adamson, Gill, and Pople */

static void
gga_xc_edf1_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 1;
  p->mix->gga_n = 3;
  XC(mix_func_alloc)(p->mix);

  XC(lda_x_init)(&p->mix->lda_mix[0], p->nspin, 3, XC_NON_RELATIVISTIC);
  p->mix->lda_coef[0] = (1.030952 - 10.4017 + 8.44793);

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_B88, p->nspin);
  XC(gga_x_b88_set_params)(&p->mix->gga_mix[0], 0.0035);
  p->mix->gga_coef[0] = 10.4017;

  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_X_B88, p->nspin);
  XC(gga_x_b88_set_params)(&p->mix->gga_mix[1], 0.0042);
  p->mix->gga_coef[1] = -8.44793;

  XC(gga_init)(&p->mix->gga_mix[2], XC_GGA_C_LYP, p->nspin);
  XC(gga_c_lyp_set_params)(&p->mix->gga_mix[2], 0.055, 0.158, 0.25, 0.3505);
  p->mix->gga_coef[2] = 1.0;
}

const XC(func_info_type) XC(func_info_gga_xc_edf1) = {
  XC_GGA_XC_EDF1,
  XC_EXCHANGE_CORRELATION,
  "EDF1",
  XC_FAMILY_GGA,
  "RD Adamson, PMW Gill, and JA Pople, Chem. Phys. Lett. 284 6 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_edf1_init, 
  NULL, NULL, NULL
};
