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

#define XC_HYB_GGA_XC_B3PW91 401 /* The original hybrid proposed by Becke */

void gga_xc_b3pw91_init(void *p_)
{
  const FLOAT a0 = 0.20, ax = 0.72, ac = 0.81;

  xc_hyb_gga_type *p = (xc_hyb_gga_type *)p_;

  p->lda_n = 2;
  p->gga_n = 2;

  xc_hyb_gga_alloc(p);

  p->exx_coef = a0;

  xc_lda_x_init(p->lda_aux[0], p->nspin, 3, XC_NON_RELATIVISTIC);
  p->lda_coef[0] = 1.0 - a0 - ax;
  xc_lda_init  (p->lda_aux[1], XC_LDA_C_PW, p->nspin);
  p->lda_coef[1] = 1.0 - ac;

  xc_gga_init(p->gga_aux[0], XC_GGA_X_B88, p->nspin);
  p->gga_coef[0] = ax;
  xc_gga_init(p->gga_aux[1], XC_GGA_C_PW91, p->nspin);
  p->gga_coef[1] = ac;
}


const xc_func_info_type func_info_hyb_gga_xc_b3pw91 = {
  XC_HYB_GGA_XC_B3PW91,
  XC_EXCHANGE_CORRELATION,
  "B3PW91",
  XC_FAMILY_HYB_GGA,
  "AD Becke, J. Chem. Phys. 98, 5648 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3pw91_init,
  NULL, 
  NULL,
  NULL /* this is taken care by the generic routine */
};
