/*
 Copyright (C) 2014 Susi Lehtola

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

#define XC_HYB_MGGA_XC_TPSSH       457 /*    TPSS hybrid */
#define XC_HYB_MGGA_XC_REVTPSSH    458 /* revTPSS hybrid */

static void
hyb_mgga_xc_tpssh_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_TPSS, XC_MGGA_C_TPSS};
  static FLOAT funcs_coef[2] = {0.9, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.10;
}

void 
XC(hyb_mgga_xc_tpssh_set_params)(XC(func_type) *p, FLOAT alpha)
{
  assert(alpha>=0 && alpha<=1.0);

  p->cam_alpha   = alpha;
  p->mix_coef[0] = 1.0 - alpha;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_tpssh) = {
  XC_HYB_MGGA_XC_TPSSH,
  XC_EXCHANGE_CORRELATION,
  "TPSSh",
  XC_FAMILY_HYB_MGGA,
  "V. N. Staroverov, G. E. Scuseria, J. Tao and J. P. Perdew, J. Chem. Phys. 119, 12129 (2003)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_mgga_xc_tpssh_init,
  NULL, NULL, NULL, NULL /* this is taken care by the generic routine */
};


static void
hyb_mgga_xc_revtpssh_init(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_REVTPSS, XC_MGGA_C_REVTPSS};
  static FLOAT funcs_coef[2] = {0.9, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.10;
}

void 
XC(hyb_mgga_xc_revtpssh_set_params)(XC(func_type) *p, FLOAT alpha)
{
  assert(alpha>=0 && alpha<=1.0);

  p->cam_alpha   = alpha;
  p->mix_coef[0] = 1.0 - alpha;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_revtpssh) = {
  XC_HYB_MGGA_XC_REVTPSSH,
  XC_EXCHANGE_CORRELATION,
  "revTPSSh",
  XC_FAMILY_HYB_MGGA,
  "G. I. Csonka, J. P. Perdew, and A. Ruzsinszky, J. Chem. Theory Comput. 6, 3688 (2010)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_mgga_xc_revtpssh_init,
  NULL, NULL, NULL, NULL /* this is taken care by the generic routine */
};
