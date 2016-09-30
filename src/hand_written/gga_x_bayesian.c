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

void XC(gga_x_bayesian_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT ss, aux, f0, f02, df0, d2f0, d3f0;

  ss = X2S*x;

  aux = 1.0 + ss;
  f0  = ss/aux;
  f02 = f0*f0;

  *f = theta0 + f02*(theta1 + f02*theta2);

  if(order < 1) return;

  df0 = 1.0/(aux*aux);

  *dfdx  = 2.0*f0*(theta1 + 2.0*theta2*f02)*df0;
  *dfdx *= X2S;

  if(order < 2) return;

  d2f0 = -2.0*df0/aux;

  *d2fdx2  = 2.0*(theta1 + 6.0*theta2*f02)*df0*df0 + 2.0*f0*(theta1 + 2.0*theta2*f02)*d2f0;
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3f0 = -3.0*d2f0/aux;

  *d3fdx3  = 24.0*theta2*f0*df0*df0*df0 + 6.0*(theta1 + 6.0*theta2*f02)*df0*d2f0 +
    2.0*f0*(theta1 + 2.0*theta2*f02)*d3f0;
  *d3fdx3 *= X2S*X2S*X2S;
  
}

#define func XC(gga_x_bayesian_enhance)
