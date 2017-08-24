/*
 Copyright (C) 2006-2009 M.A.L. Marques

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

#define XC_LDA_C_1D_LOOS          26 /* P-F Loos correlation LDA     */

#include "maple2c/lda_c_1d_loos.c"

#define func maple2c_func
#define XC_DIMENSIONS 1
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_1d_loos) = {
  XC_LDA_C_1D_LOOS,
  XC_CORRELATION,
  "P-F Loos correlation LDA",
  XC_FAMILY_LDA,
  {&xc_ref_Loos2013_064108, NULL, NULL, NULL, NULL},
  XC_FLAGS_1D |  XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0,
  0, NULL, NULL,
  NULL, NULL,
  work_lda, NULL, NULL
};
