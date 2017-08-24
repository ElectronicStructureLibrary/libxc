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

#define XC_MGGA_X_SCAN          263 /* SCAN exchange of Sun, Ruzsinszky, and Perdew  */
#define XC_HYB_MGGA_X_SCAN0     264 /* SCAN hybrid exchange */

/* WARNING THIS CODE WILL BE DELETED */
void
XC(mgga_x_scan_falpha)(int order, FLOAT a, FLOAT c1, FLOAT c2, FLOAT dd, FLOAT *f, FLOAT *dfda)
{
  /* exponentials are truncated */
  const FLOAT logeps =  LOG(FLOAT_EPSILON);
  FLOAT thr1, thr2;
  FLOAT c1exp, c2exp, ooma;

  thr1  = -logeps/(c1 - logeps);
  thr2  = 1.0 - c2/logeps;

  ooma = 1.0/(1.0 - a);

  c1exp = (a >= thr1) ? 0.0 : EXP(-c1*a*ooma);
  c2exp = (a <= thr2) ? 0.0 : EXP(c2*ooma);

  *f = c1exp - dd*c2exp;

  if(order < 1) return;

  *dfda = -(c1*c1exp + dd*c2*c2exp)*ooma*ooma;
}

#include "maple2c/mgga_x_scan.c"

#define func maple2c_func
#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_scan) = {
  XC_MGGA_X_SCAN,
  XC_EXCHANGE,
  "SCAN exchange of Sun, Ruzsinszky, and Perdew",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  NULL, NULL, NULL, NULL,
  work_mgga_x,
};

static void
hyb_mgga_x_scan0_init(XC(func_type) *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_SCAN};
  static FLOAT funcs_coef[1] = {1.0 - 0.25};

  XC(mix_init)(p, 1, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}


const XC(func_info_type) XC(func_info_hyb_mgga_x_scan0) = {
  XC_HYB_MGGA_X_SCAN0,
  XC_EXCHANGE,
  "SCAN hybrid exchange (SCAN0)",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Hui2016_044114, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32,
  0, NULL, NULL,
  hyb_mgga_x_scan0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
