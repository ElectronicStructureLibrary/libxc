/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_WC         118 /* Wu & Cohen */

static FLOAT wc_mu, wc_c;

static void
gga_x_wc_init(void *p_)
{
  wc_mu  = 0.2195149727645171;
  wc_c   = (146.0/2025.0)*(4.0/9.0) - (73.0/405.0)*(2.0/3.0) + (wc_mu - 10.0/81.0);
}

static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  const FLOAT kappa = 0.8040;

  FLOAT s, s2;
  FLOAT aux1, aux2, f0, df0, d2f0, dd;

  s  = X2S*x;
  s2 = s*s;
  
  aux1 = wc_mu - 10.0/81.0;
  aux2 = exp(-s2);

  f0 = kappa + 10.0/81.0*s2 + s2*aux1*aux2 + log(1.0 + wc_c*s2*s2);
  *f = 1.0 + kappa*(1.0 - kappa/f0);

  if(dfdx==NULL && d2fdx2==NULL) return; /* nothing else to do */

  df0 = 20.0/81.0*s + 2.0*s*aux1*aux2*(1.0 - s2) + 4.0*wc_c*s*s2/(1.0 + wc_c*s2*s2);

  if(dfdx!=NULL){
    *dfdx  = X2S*kappa*kappa*df0/(f0*f0);
    *ldfdx = X2S*X2S*wc_mu;
  }

  if(d2fdx2==NULL) return; /* nothing else to do */

  dd   = 1.0 + wc_c*s2*s2;
  d2f0 = 20.0/81.0 + 2.0*aux1*aux2*(1.0 - 5.0*s2 + 2.0*s2*s2)
    - 4.0*wc_c*s2*(dd - 4.0)/(dd*dd);

  *d2fdx2 = X2S*X2S*kappa*kappa/(f0*f0)*(d2f0 - 2.0*df0*df0/f0);
}


#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_wc) = {
  XC_GGA_X_WC,
  XC_EXCHANGE,
  "Wu & Cohen",
  XC_FAMILY_GGA,
  "Z Wu and RE Cohen, Phys. Rev. B 73, 235116 (2006)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  gga_x_wc_init, 
  NULL, NULL,
  work_gga_x
};
