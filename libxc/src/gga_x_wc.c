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

/* 
     the variable used in the original PBE paper is s = |grad n|/(2 k_f n)
     while here we use x = |grad n_s|/n_s^(4/3). Therefore, the value of s
     we use here is the original x multiplied by 

       1/2^(1/3) * 1/(2*(3 pi^2)^(1/3))

     where the first term comes from using the spin densities, and the second
     from the constants of k_f = (3 pi^2 n)^(1/3)
*/

static FLOAT wc_mu;
static FLOAT wc_c, wc_fac;

static void
gga_x_wc_init(void *p_)
{
  wc_mu  = 0.2195149727645171;
  wc_fac = 1.0/(2.0*POW(6.0*M_PI*M_PI, 1.0/3.0));
  wc_c   = (146.0/2025.0)*(4.0/9.0) - (73.0/405.0)*(2.0/3.0) + (wc_mu - 10.0/81.0);
}

static inline void 
func(XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx)
{
  const FLOAT kappa = 0.8040;

  FLOAT s, s2, xx, dxx;
  FLOAT aux1, aux2, dd;

  s  = wc_fac*x;
  s2 = s*s;
  
  aux1 = wc_mu - 10.0/81.0;
  aux2 = exp(-s2);

  xx = 10.0/81.0*s2 + s2*aux1*aux2 + log(1.0 + wc_c*s2*s2);
  dd = 1.0/(kappa + xx);

  *f = 1.0 + kappa*(1.0 - kappa*dd);

  dxx  = 20.0/81.0*s + 2.0*s*aux1*aux2*(1.0 - s2) + 4.0*wc_c*s*s2/(1.0 + wc_c*s2*s2);
  dxx *= wc_fac; /* convert d/ds in d/dx */

  *dfdx  = dxx * kappa*kappa * dd*dd;
  *ldfdx = 10.0/81.0;
}


#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_wc) = {
  XC_GGA_X_WC,
  XC_EXCHANGE,
  "Wu & Cohen",
  XC_FAMILY_GGA,
  "Z Wu and RE Cohen, Phys. Rev. B 73, 235116 (2006)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_x_wc_init, 
  NULL, NULL,
  work_gga_x
};
