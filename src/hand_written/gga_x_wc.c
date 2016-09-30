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

void XC(gga_x_wc_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT ss, ss2;
  FLOAT aux1, aux2, f0, df0, d2f0, d3f0, dd;

  ss  = X2S*x;
  ss2 = ss*ss;
  
  aux1 = mu - MU_GE;
  aux2 = EXP(-ss2);

  f0 = kappa + MU_GE*ss2 + ss2*aux1*aux2 + LOG(1.0 + c*ss2*ss2);
  *f = 1.0 + kappa*(1.0 - kappa/f0);

  if(order < 1) return;

  dd   = 1.0 + c*ss2*ss2;
  df0 = 2.0*MU_GE*ss + 2.0*ss*aux1*aux2*(1.0 - ss2) + 4.0*c*ss*ss2/dd;

  *dfdx  = kappa*kappa*df0/(f0*f0);
  *dfdx *= X2S;

  if(order < 2) return;

  d2f0 = 2.0*MU_GE + 2.0*aux1*aux2*(1.0 - 5.0*ss2 + 2.0*ss2*ss2)
    - 4.0*c*ss2*(dd - 4.0)/(dd*dd);

  *d2fdx2  = -kappa*kappa*(2.0*df0*df0 - d2f0*f0)/(f0*f0*f0);
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3f0 = -4.0*aux1*aux2*ss*(6.0 - 9.0*ss2 + 2.0*ss2*ss2) +
    8.0*c*ss*(3.0 + c*ss2*ss2*(c*ss2*ss2 - 12.0))/(dd*dd*dd);

  *d3fdx3  = kappa*kappa*(6.0*df0*df0*df0 - 6.0*f0*df0*d2f0 + f0*f0*d3f0)/(f0*f0*f0*f0);
  *d3fdx3 *= X2S*X2S*X2S;
 
}

#define func XC(gga_x_wc_enhance)
