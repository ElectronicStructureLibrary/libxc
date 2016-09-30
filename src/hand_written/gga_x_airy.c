/*
 Copyright (C) 2008 Georg Madsen

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

void XC(gga_x_airy_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT ss, ssa2, ssa6, ssa8, ssa10;
  FLOAT   den1,   aux1,   num1,   den2,   num2;
  FLOAT  dden1,  daux1,  dnum1,  dden2,  dnum2;
  FLOAT d2den1, d2aux1, d2num1, d2den2, d2num2;
  FLOAT d3den1, d3aux1, d3num1, d3den2, d3num2;

  ss = X2S*x;
  ssa2  = POW(ss, a2);
  ssa6  = POW(ss, a6);
  ssa8  = POW(ss, a8);
  ssa10 = POW(ss, a10);

  num1 = a1*ssa2;
  aux1 = 1.0 + a3*ssa2;
  den1 = POW(aux1, a4);

  *f = num1/den1;

  if(p->info->number == XC_GGA_X_AIRY){
    num2 = 1.0 - a5*ssa6 + a7*ssa8;
    den2 = 1.0 + a9*ssa10;

    *f += num2/den2;
  }

  if(order < 1) return;

  dnum1 = a2*num1/ss;
  daux1 = a3*a2*ssa2/ss;
  dden1 = a4*daux1*den1/aux1;

  *dfdx  = DFRACTION(num1, dnum1, den1, dden1);

  if(p->info->number == XC_GGA_X_AIRY){
    dnum2 = -a5*a6*ssa6/ss + a7*a8*ssa8/ss;
    dden2 = a9*a10*ssa10/ss;

    *dfdx += DFRACTION(num2, dnum2, den2, dden2);
  }

  *dfdx *= X2S;

  if(order < 2) return;

  d2num1 = (a2 - 1.0)*dnum1/ss;
  d2aux1 = (a2 - 1.0)*daux1/ss;
  d2den1 = a4*den1/(aux1*aux1) * ((a4 - 1.0)*daux1*daux1 + aux1*d2aux1);

  *d2fdx2 = D2FRACTION(num1, dnum1, d2num1, den1, dden1, d2den1);

  if(p->info->number == XC_GGA_X_AIRY){
    d2num2 = -a5*a6*(a6 - 1.0)*ssa6/(ss*ss) + a7*a8*(a8 - 1.0)*ssa8/(ss*ss);
    d2den2 = (a10 - 1.0)*dden2/ss;

    *d2fdx2 += D2FRACTION(num2, dnum2, d2num2, den2, dden2, d2den2);
  }
  
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3num1 = (a2 - 2.0)*d2num1/ss;
  d3aux1 = (a2 - 1.0)*d2aux1/ss;
  d3den1 = a4*den1/(aux1*aux1*aux1) * ((a4 - 1.0)*daux1*((a4 - 2.0)*daux1*daux1 + 3.0*aux1*d2aux1) + aux1*aux1*d3aux1);

  *d3fdx3 = D3FRACTION(num1, dnum1, d2num1, d3num1, den1, dden1, d2den1, d3den1);

  if(p->info->number == XC_GGA_X_AIRY){
    d3num2 = -a5*a6*(a6 - 1.0)*(a6 - 2.0)*ssa6/(ss*ss*ss) + a7*a8*(a8 - 1.0)*(a8 - 2.0)*ssa8/(ss*ss*ss);
    d3den2 = (a10 - 2.0)*d2den2/ss;

    *d3fdx3 += D3FRACTION(num2, dnum2, d2num2, d3num2, den2, dden2, d2den2, d3den2);
  }

  *d3fdx3 *= X2S*X2S*X2S;
}

#define func XC(gga_x_airy_enhance)
