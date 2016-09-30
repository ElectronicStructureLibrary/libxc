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

void XC(gga_x_rpbe_enhance) 
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT kappa, mu, f0, df0, d2f0, d3f0;

  assert(p->params != NULL);
  kappa = ((gga_x_rpbe_params *) (p->params))->kappa;
  mu    = ((gga_x_rpbe_params *) (p->params))->mu*X2S*X2S;

  f0 = EXP(-mu*x*x/kappa);
  *f = 1.0 + kappa*(1.0 - f0);

  if(order < 1) return;

  df0 = -2.0*x*mu/kappa*f0;
  
  *dfdx  = -kappa*df0;

  if(order < 2) return;

  d2f0    = -2.0*mu*f0*(kappa - 2.0*x*x*mu)/(kappa*kappa);
  *d2fdx2 = -kappa*d2f0;

  if(order < 3) return;

  d3f0    = 4.0*mu*mu*f0*x*(3.0*kappa - 2.0*mu*x*x)/(kappa*kappa*kappa);
  *d3fdx3 = -kappa*d3f0;
}

#define func XC(gga_x_rpbe_enhance)
