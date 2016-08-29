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

#define XC_GGA_X_SG4         533 /* Semiclassical GGA at fourth order */


void XC(gga_x_sg4_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static FLOAT mu1 = 0.042, mu2 = 0.218; /* mu1 + mu2 = 0.26 */
  static FLOAT k2 = 0.24371282051282051282; /* k2 = -mu2^2/nu_MGE4, nu_MGE4 = -0.195 */
  static FLOAT k1 = 0.56028717948717948718; /* k1 + k2 = 0.804 */

  FLOAT ss, ssa, ssa2, ssb, dena, ddena, d2dena, d3dena;

  ss  = X2S*x;

  ssa  = mu1*ss*ss/k1;
  ssa2 = ssa*ssa;
  dena = 1.0 - ssa2*ssa2*ssa;

  ssb = mu2*ss*ss/k2;

  *f = 1.0 + k1 + k2 - k1*(1.0 - ssa)/dena - k2/(1.0 + ssb);

  if(order < 1) return;

  ddena = -10.0*ssa2*ssa2*ss*mu1/k1;

  *dfdx  = 
    - k1*DFRACTION(1.0 - ssa, -2.0*ss*mu1/k1, dena, ddena)
    - k2*DFRACTION(1.0, 0.0, 1.0 + ssb, 2.0*ss*mu2/k2);

  *dfdx *= X2S;

  if(order < 2) return;

  d2dena = -90.0*ssa2*ssa2*mu1/k1;

  *d2fdx2 = 
    - k1*D2FRACTION(1.0 - ssa, -2.0*ss*mu1/k1, -2.0*mu1/k1, dena, ddena, d2dena)
    - k2*D2FRACTION(1.0, 0.0, 0.0, 1.0 + ssb, 2.0*ss*mu2/k2, 2.0*mu2/k2);

  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3dena = -90.0*8.0*ssa2*ssa*ss*mu1*mu1/(k1*k1);

  *d3fdx3 =     
    - k1*D3FRACTION(1.0 - ssa, -2.0*ss*mu1/k1, -2.0*mu1/k1, 0.0, dena, ddena, d2dena, d3dena)
    - k2*D3FRACTION(1.0, 0.0, 0.0, 0.0, 1.0 + ssb, 2.0*ss*mu2/k2, 2.0*mu2/k2, 0.0);

  *d3fdx3 *= X2S*X2S*X2S;
}


#define func XC(gga_x_sg4_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_sg4) = {
  XC_GGA_X_SG4,
  XC_EXCHANGE,
  "Semiclassical GGA at fourth order",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2016_045126, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};
