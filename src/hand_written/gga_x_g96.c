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

void
XC(gga_x_g96_enhance)(const XC(func_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static const FLOAT c1 = 1.0/137.0;
  FLOAT sx = SQRT(x);

  *f     = 1.0 + c1/X_FACTOR_C*x*sx;

  if(order < 1) return;

  *dfdx  = 3.0*c1/(2.0*X_FACTOR_C)*sx;

  if(order < 2) return;

  *d2fdx2 = 3.0*c1/(4.0*X_FACTOR_C*sx);

  if(order < 2) return;

  *d3fdx3 = -3.0*c1/(8.0*X_FACTOR_C*x*sx);
}

#define func XC(gga_x_g96_enhance)
