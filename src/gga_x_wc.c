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
#include <assert.h>
#include "util.h"

#define XC_GGA_X_WC         118 /* Wu & Cohen */

static FLOAT wc_mu, wc_c;

static void
gga_x_wc_init(XC(func_type) *p_)
{
  wc_mu  = 0.2195149727645171;
  wc_c   = (146.0/2025.0)*(4.0/9.0) - (73.0/405.0)*(2.0/3.0) + (wc_mu - 10.0/81.0);
}

void XC(gga_x_wc_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  const FLOAT kappa = 0.8040;

  FLOAT ss, ss2;
  FLOAT aux1, aux2, f0, df0, d2f0, d3f0, dd;

  ss  = X2S*x;
  ss2 = ss*ss;
  
  aux1 = wc_mu - 10.0/81.0;
  aux2 = EXP(-ss2);

  f0 = kappa + 10.0/81.0*ss2 + ss2*aux1*aux2 + LOG(1.0 + wc_c*ss2*ss2);
  *f = 1.0 + kappa*(1.0 - kappa/f0);

  if(order < 1) return;

  dd   = 1.0 + wc_c*ss2*ss2;
  df0 = 20.0/81.0*ss + 2.0*ss*aux1*aux2*(1.0 - ss2) + 4.0*wc_c*ss*ss2/dd;

  *dfdx  = kappa*kappa*df0/(f0*f0);
  *dfdx *= X2S;

  if(order < 2) return;

  d2f0 = 20.0/81.0 + 2.0*aux1*aux2*(1.0 - 5.0*ss2 + 2.0*ss2*ss2)
    - 4.0*wc_c*ss2*(dd - 4.0)/(dd*dd);

  *d2fdx2  = -kappa*kappa*(2.0*df0*df0 - d2f0*f0)/(f0*f0*f0);
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3f0 = -4.0*aux1*aux2*ss*(6.0 - 9.0*ss2 + 2.0*ss2*ss2) +
    8.0*wc_c*ss*(3.0 + wc_c*ss2*ss2*(wc_c*ss2*ss2 - 12.0))/(dd*dd*dd);

  *d3fdx3  = kappa*kappa*(6.0*df0*df0*df0 - 6.0*f0*df0*d2f0 + f0*f0*d3f0)/(f0*f0*f0*f0);
  *d3fdx3 *= X2S*X2S*X2S;
 
}

#define func XC(gga_x_wc_enhance)
#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_wc) = {
  XC_GGA_X_WC,
  XC_EXCHANGE,
  "Wu & Cohen",
  XC_FAMILY_GGA,
  "Z Wu and RE Cohen, Phys. Rev. B 73, 235116 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_wc_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};
