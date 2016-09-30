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
