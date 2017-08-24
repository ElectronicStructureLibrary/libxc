/*
 Copyright (C) 2017 M.A.L. Marques

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

#define XC_LDA_X_ERF   546   /* Attenuated exchange LDA (erf) */

static void lda_x_erf_init(XC(func_type) *p)
{
  /* initialize omega to something reasonable */
  p->cam_omega = 0.3;
}


#include "maple2c/lda_x_erf.c"

#define func maple2c_func
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_x_erf) = {
  XC_LDA_X_ERF,
  XC_EXCHANGE,
  "Attenuated exchange LDA (erf)",
  XC_FAMILY_LDA,
  {&xc_ref_Toulouse2004_1047, &xc_ref_Tawada2004_8425, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-29, 0.0, 0.0,
  0, NULL, NULL,
  lda_x_erf_init, NULL, 
  work_lda, NULL, NULL
};
