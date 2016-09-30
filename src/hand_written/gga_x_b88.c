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
XC(gga_x_b88_enhance)(const XC(func_type) *p, int order, FLOAT x, 
		      FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT x2, aux1, aux2, f1, f2, df1, df2, d2f1, d2f2, d3f1, d3f2, dd;
  gga_x_b88_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b88_params *) (p->params);

  x2 = x*x;

  f1 = params->beta/X_FACTOR_C*x2;
  f2 = 1.0 + params->gamma*params->beta*x*ASINH(x);
  *f = 1.0 + f1/f2;

  if(p->info->number == XC_GGA_K_THAKKAR){
    dd  = 1.0/(1.0 + 2.0*CBRT(4.0)*x);
    *f += -0.072*x*dd;
  }

  if(order < 1) return;

  aux1 = 1.0 + x2;
  aux2 = SQRT(aux1);

  df1 = 2.0*params->beta/X_FACTOR_C*x;
  df2 = params->gamma*params->beta*(ASINH(x) + x/aux2);

  *dfdx = (df1*f2 - f1*df2)/(f2*f2);

  if(p->info->number == XC_GGA_K_THAKKAR)
    *dfdx += -0.072*dd*dd;
    
  if(order < 2) return;

  d2f1 = 2.0*params->beta/X_FACTOR_C;
  d2f2 = params->gamma*params->beta*(2.0 + x2)/(aux1*aux2);

  *d2fdx2 = (2.0*f1*df2*df2 + d2f1*f2*f2 - f2*(2.0*df1*df2 + f1*d2f2))/(f2*f2*f2);

  if(p->info->number == XC_GGA_K_THAKKAR)
    *d2fdx2 += 0.072*4.0*CBRT(4.0)*dd*dd*dd;

  if(order < 3) return;

  d3f1 = 0.0;
  d3f2 = -params->beta*params->gamma*x*(4.0 + x2)/(aux1*aux1*aux2);

  *d3fdx3 = (-6.0*f1*df2*df2*df2 + 6.0*f2*df2*(df1*df2 + f1*d2f2) + f2*f2*f2*d3f1 - f2*f2*(3.0*df2*d2f1 + 3.0*df1*d2f2 + f1*d3f2))/(f2*f2*f2*f2);

  if(p->info->number == XC_GGA_K_THAKKAR)
    *d3fdx3 += -0.072*24.0*CBRT(4.0)*CBRT(4.0)*dd*dd*dd*dd;
}

#define func XC(gga_x_b88_enhance)
