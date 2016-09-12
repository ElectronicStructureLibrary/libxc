/*
 Copyright (C) 2015 M.A.L. Marques, Markus Patzold

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
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_GG99  535 /* Gilbert and Gill 1999 */

inline static void 
r_x(int order, FLOAT x, FLOAT *r, FLOAT *dr, FLOAT *d2r, FLOAT *d3r,FLOAT *woy,FLOAT *woz)
{
  static const FLOAT
    a1 = 4.0*M_SQRT3*M_PI*M_PI*M_PI;

  FLOAT a2, x2, x4, x6, aux1, aux2, daux1, daux2, num, den, dd, dnum, dden;

  a2 = SQRT(3.0/(2.0*a1));

  x2 = x*x;
  x4 = x2*x2;
  x6 = x2*x4;

  aux1 = a1 + SQRT(a1*a1 - x6);
  aux2 = CBRT(aux1);

  num = x*a2*SQRT(x2 + aux2*aux2);
  den = SQRT(aux2);
  
  *woy = num;
  *woz = den;

  *r = asinh(num/den);

  if(order < 1) return;

  daux1 = -3.0*x*x4/SQRT(a1*a1 - x6);
  daux2 = daux1*aux2/(3.0*aux1);

  dnum = a2*(2.0*x2 + aux2*aux2 + x*aux2*daux2)/SQRT(x2 + aux2*aux2);
  dden = daux2*den/(2.0*aux2);

  dd  = DFRACTION(num, dnum, den, dden);
  *dr = dd/SQRT(1 + num*num/(den*den));
}


void XC(gga_x_gg99_enhance)
     (const XC(func_type) *p, int order, FLOAT x, 
      FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT r, dr, d2r, d3r;
  FLOAT aux1, aux2, aux3, aux4, aux5, daux1, daux2, daux4, daux5;
  FLOAT num, den, dnum, dden, df;
  FLOAT x2,x5,x6;
  FLOAT woy,woz,ablw,ablr;
  FLOAT ablr1,ablr2,abl1,abl2;
  
  r_x(order, x, &r, &dr, &d2r, &d3r, &woy, &woz);
  
  aux1 = EXP(-2.0*r);

  aux2 = LOG(1.0 + aux1);
  aux3 = 1.0/cosh(r);
  aux4 = POW(aux3, 2.0/3.0);
  aux5 = XC(dilogarithm)(-aux1);

  num = -M_PI*M_PI + 12.0*r*aux2 - 12.0*aux5;
  den = 2.0*M_CBRT3*M_PI*r*aux4;

  *f = num/den;
    
  if(order < 1) return;

  daux1 = -2.0*aux1;
  daux2 = daux1/(1.0 + aux1);
  daux4 = -2.0/3.0*aux4*tanh(r);
  daux5 = -aux2*daux1/aux1;

  dnum = 12.0*(aux2 + r*daux2) - 12.0*daux5;
  dden = 2.0*M_CBRT3*M_PI*(aux4 + r*daux4);

  *dfdx = DFRACTION(num, dnum, den, dden)*dr;
}

#define func XC(gga_x_gg99_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_gg99) = {
  XC_GGA_X_GG99,
  XC_EXCHANGE,
  "Gilbert and Gill 1999",
  XC_FAMILY_GGA,
  {&xc_ref_Gilbert1999_511, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, work_gga_x, NULL
};
