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

#define XC_GGA_X_C09X         158 /* C09x to be used with the VdW of Rutgers-Chalmers     */

void XC(gga_x_c09x_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static FLOAT mu = 0.0617, kappa = 1.245, alpha = 0.0483;

  FLOAT ss, ss2, aux, aux2;

  ss  = X2S*x;
  ss2 = ss*ss;

  aux  = EXP(-0.5*alpha*ss2);
  aux2 = aux*aux;

  *f = 1.0 + mu*ss2*aux2 + kappa*(1.0 - aux);

  if(order < 1) return;

  *dfdx  = 2.0*ss*mu*aux2*(1.0 - alpha*ss2) + alpha*kappa*ss*aux;
  *dfdx *= X2S;

  if(order < 2) return;

  *d2fdx2  = 2.0*mu*aux2*(1.0 + alpha*ss2*(2.0*alpha*ss2 - 5.0)) + alpha*kappa*aux*(1.0 - alpha*ss2);
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;
  
  *d3fdx3  = alpha*ss*(-4.0*mu*aux2*(6.0 + alpha*ss2*(2.0*alpha*ss2 - 9.0)) + alpha*kappa*aux*(alpha*ss2 - 3.0));
  *d3fdx3 *= X2S*X2S*X2S;
  
}


#define func XC(gga_x_c09x_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_c09x) = {
  XC_GGA_X_C09X,
  XC_EXCHANGE,
  "C09x to be used with the VdW of Rutgers-Chalmers",
  XC_FAMILY_GGA,
  "VR Cooper, Phys. Rev. B 81, 161104(R) (2010)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};
