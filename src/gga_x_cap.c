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

#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_CAP         270 /* Correct Asymptotic Potential */
#define XC_HYB_GGA_XC_CAP0   477 /* Correct Asymptotic Potential hybrid */

void XC(gga_x_cap_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  const FLOAT Ax=-3.0/4.0*CBRT(3.0/M_PI);
  const FLOAT mu=0.2195149727645171;
  const FLOAT alpha=-Ax*mu;
  const FLOAT c=alpha/CBRT(3*M_PI*M_PI);
  const FLOAT alphaoAx=-mu;
  
  FLOAT s;
  FLOAT ops;
  FLOAT logops;
  FLOAT opclogops;

  s  = X2S*x;

  ops = 1+s;
  logops = LOG(ops);
  opclogops = 1 + c*logops;

  /* eqn (14) */
  *f = 1.0 - alphaoAx * (s*logops) / opclogops;

  if(order < 1) return;

  /* eqn (A3) */
  *dfdx = -alphaoAx * (s + ops*logops*opclogops) / (ops*opclogops*opclogops) * X2S;
  
  if(order < 2) return;

  /* eqn (A4) */
  *d2fdx2 = -alphaoAx * (1 + (1.0-2.0*c)*s + c*(2+s)*logops ) / (ops*ops*logops*logops*logops) * X2S * X2S;
}

#define func XC(gga_x_cap_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_cap) = {
  XC_GGA_X_CAP,
  XC_EXCHANGE,
  "Correct Asymptotic Potential",
  XC_FAMILY_GGA,
  {&xc_ref_Carmona2015_054105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};

void
XC(hyb_gga_xc_cap0_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_CAP, XC_GGA_C_PBE};
  static FLOAT funcs_coef[2] = {0.8, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  /* C functional is PBE C with β = (3/4)β PBE */
  XC(gga_c_pbe_set_params)(p->func_aux[1],0.75*0.06672455060314922);
  
  p->cam_alpha = 0.20;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_cap0) = {
  XC_HYB_GGA_XC_CAP0,
  XC_EXCHANGE_CORRELATION,
  "Correct Asymptotic Potential hybrid",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Carmona2016_120, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  XC(hyb_gga_xc_cap0_init),
  NULL, NULL, NULL, NULL
};
