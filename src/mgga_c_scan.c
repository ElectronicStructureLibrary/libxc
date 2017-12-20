/*
 Copyright (C) 2016 Susi Lehtola

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

#define XC_MGGA_C_SCAN          267 /* SCAN correlation */
#define XC_MGGA_C_SCAN_RVV10    292 /* SCAN correlation + rVV10 correlation */

#include "maple2c/mgga_c_scan.c"

#define func maple2c_func
#include "work_mgga_c.c"

const xc_func_info_type xc_func_info_mgga_c_scan = {
  XC_MGGA_C_SCAN,
  XC_CORRELATION,
  "SCAN correlation of Sun, Ruzsinszky, and Perdew",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-26,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, NULL, work_mgga_c,
};


static void
mgga_c_scan_rvv10_init(xc_func_type *p)
{
  static int   funcs_id  [1] = {XC_MGGA_C_SCAN};
  static double funcs_coef[1] = {1.0};

  xc_mix_init(p, 1, funcs_id, funcs_coef);

  p->nlc_b = 15.7;
  p->nlc_C = 0.0093;
}

const xc_func_info_type xc_func_info_mgga_c_scan_rvv10 = {
  XC_MGGA_C_SCAN_RVV10,
  XC_CORRELATION,
  "SCAN+rVV10 correlation",
  XC_FAMILY_MGGA,
  {&xc_ref_Peng2016_041005, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_VV10 | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32,
  0, NULL, NULL,
  mgga_c_scan_rvv10_init,
  NULL, NULL, NULL, NULL
};
