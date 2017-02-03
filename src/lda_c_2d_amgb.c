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


#include "util.h"

/************************************************************************
 Correlation energy per particle and potentials for a homogeneous electron
 gas in 2D, as parametrized by Attaccalite et al.
************************************************************************/

#define XC_LDA_C_2D_AMGB  15   /* Attaccalite et al             */

#include "maple2c/lda_c_2d_amgb.c"

#define func maple2c_func
#define XC_DIMENSIONS 2
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_2d_amgb) = {
  XC_LDA_C_2D_AMGB,
  XC_CORRELATION,
  "AMGB (for 2D systems)",
  XC_FAMILY_LDA,
  {&xc_ref_Attaccalite2002_256601, NULL, NULL, NULL, NULL},
  XC_FLAGS_2D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-9, 0.0, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, NULL,
  work_lda, NULL, NULL
};
