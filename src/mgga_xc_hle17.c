/*
 Copyright (C) 2017 Susi Lehtola

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

#define XC_MGGA_XC_HLE17      288  /* high local exchange 2017   */

static void
mgga_xc_hle17_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
  static FLOAT funcs_coef[2] = {1.25, 0.5};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
}

const XC(func_info_type) XC(func_info_mgga_xc_hle17) = {
  XC_MGGA_XC_HLE17,
  XC_EXCHANGE_CORRELATION,
  "high local exchange 2017",
  XC_FAMILY_MGGA,
  {&xc_ref_Verma2017, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  mgga_xc_hle17_init,
  NULL, NULL, NULL, NULL
};
