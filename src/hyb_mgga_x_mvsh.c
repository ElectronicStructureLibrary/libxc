/*
 Copyright (C) 2015 Susi Lehtola

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

#define XC_HYB_MGGA_X_MVSH     474 /* MVSh hybrid */

static void
hyb_mgga_x_mvsh_init(XC(func_type) *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_MVS};
  static FLOAT funcs_coef[1] = {0.75};

  XC(mix_init)(p, 1, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}

void 
XC(hyb_mgga_x_mvsh_set_params)(XC(func_type) *p, FLOAT alpha)
{
  assert(alpha>=0 && alpha<=1.0);

  p->cam_alpha   = alpha;
  p->mix_coef[0] = 1.0 - alpha;
}

const XC(func_info_type) XC(func_info_hyb_mgga_x_mvsh) = {
  XC_HYB_MGGA_X_MVSH,
  XC_EXCHANGE,
  "MVSh hybrid exchange functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Sun2015_685, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_mgga_x_mvsh_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
