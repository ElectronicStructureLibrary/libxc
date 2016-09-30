/*
 Copyright (C) 2006-2014 L. Talirz, M.A.L. Marques

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

void XC(gga_x_b86_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT dd, ddd, d2dd;
  FLOAT f1, f2, df1, df2, d2f1, d2f2, d3f2;

  gga_x_b86_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b86_params *) (p->params);

  dd = 1.0 + params->gamma*x*x;
  f1 = params->beta*x*x;
  f2 = POW(dd, params->omega);

  *f = 1.0 + f1/f2;
  
  if(order < 1) return;

  ddd = 2.0*params->gamma*x;
  df1 = 2.0*params->beta *x;
  df2 = params->omega*ddd*f2/dd;

  *dfdx  = DFRACTION(f1, df1, f2, df2);

  if(order < 2) return;

  d2dd = 2.0*params->gamma;
  d2f1 = 2.0*params->beta;
  d2f2 = params->omega*f2/(dd*dd)*(d2dd*dd + (params->omega - 1.0)*ddd*ddd);

  *d2fdx2 = D2FRACTION(f1, df1, d2f1, f2, df2, d2f2);

  if(order < 3) return;

  d3f2 = params->omega*(params->omega - 1.0)*ddd*f2/(dd*dd*dd)*(3.0*d2dd*dd + (params->omega - 2.0)*ddd*ddd);

  *d3fdx3 = D3FRACTION(f1, df1, d2f1, 0.0, f2, df2, d2f2, d3f2);
}

#define func XC(gga_x_b86_enhance)
