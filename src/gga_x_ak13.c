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
#include "util.h"

#define XC_GGA_X_AK13  56 /* Armiento & Kuemmel 2013 */

static const FLOAT B1 =  1.74959015598863046792081721182; /* 3*muGE/5 + 8 pi/15 */
static const FLOAT B2 = -1.62613336586517367779736042170; /* muGE - B1 */

FLOAT XC(gga_ak13_get_asymptotic) (FLOAT homo)
{
  FLOAT Qx, aa, aa2, factor;

  Qx = SQRT(2.0)*B1/(3.0*CBRT(3.0*M_PI*M_PI));

  aa  = X_FACTOR_C*Qx;
  aa2 = aa*aa;

  factor = (homo < 0.0) ? -1.0 : 1.0;
    
  return (aa2/2.0)*(1.0 + factor*SQRT(1.0 - 4.0*homo/aa2));
}


void XC(gga_x_ak13_enhance) 
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT ss, den, f1, f2;

  ss  = X2S*x;

  f1 = LOG(1.0 + ss);
  f2 = LOG(1.0 + f1);

  *f = 1.0 + B1*ss*f1 + B2*ss*f2;

  if(order < 1) return;

  den = (1.0 + ss)*(1.0 + f1);

  *dfdx  = B1*f1 + ss*(B1 + B2 + B1*f1)/den + B2*f2;
  *dfdx *= X2S;

  if(order < 2) return;
  
  *d2fdx2  = (2.0*B2 + B1*(2.0 + ss) + (2.0 + ss)*f1*(2.0*B1 + B2 + B1*f1))/(den*den);
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  *d3fdx3  = (B2*(ss - 6.0) - B1*(ss + 3.0) - f1*(3.0*B1*(ss + 3.0) + B2*(2.0*ss + 9.0) + (ss + 3.0)*f1*(3.0*B1 + B2 + B1*f1)))/(den*den*den);
  *d3fdx3 *= X2S*X2S*X2S;
}

#define func XC(gga_x_ak13_enhance) 
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_ak13) = {
  XC_GGA_X_AK13,
  XC_EXCHANGE,
  "Armiento & Kuemmel 2013",
  XC_FAMILY_GGA,
  "R Armiento and S Kuemmel, PRL 111, 036402 (2013)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};

