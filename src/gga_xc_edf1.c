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

#define XC_GGA_XC_EDF1 165 /* Empirical functionals from Adamson, Gill, and Pople */

static void
gga_xc_edf1_init(void *p_)
{
  static int   funcs_id  [4] = {XC_LDA_X, XC_GGA_X_B88, XC_GGA_X_B88, XC_GGA_C_LYP};
  static FLOAT funcs_coef[4] = {1.030952 - 10.4017 + 8.44793, 10.4017, -8.44793, 1.0};
  XC(gga_type) *p = (XC(gga_type) *)p_;

  gga_init_mix(p, 4, funcs_id, funcs_coef);  

  XC(gga_x_b88_set_params)(p->func_aux[1], 0.0035);
  XC(gga_x_b88_set_params)(p->func_aux[2], 0.0042);
  XC(gga_c_lyp_set_params)(p->func_aux[3], 0.055, 0.158, 0.25, 0.3505);
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
