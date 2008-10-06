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

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_PBE_SOL      116 /* Perdew, Burke & Ernzerhof exchange (solids)    */
#define XC_GGA_X_XPBE         123 /* xPBE reparametrization by Xu & Goddard         */

static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT kappa[4] = {
    0.8040,  /* original PBE */
    1.245,   /* PBE R */
    0.8040,  /* PBE sol */
    0.91954  /* xPBE */
  };

  static const FLOAT mu[4] = {
    0.2195149727645171,  /* PBE: mu = beta*pi^2/3, beta = 0.066725 */
    0.2195149727645171,  /* PBE rev: as PBE */
    10.0/81.0,           /* PBE sol */
    0.23214              /* xPBE */
  };

  FLOAT ss, f0, df0, d2f0;
  int func;

  switch(p->info->number){
  case XC_GGA_X_PBE_R:   func = 1; break;
  case XC_GGA_X_PBE_SOL: func = 2; break;
  case XC_GGA_X_XPBE:    func = 3; break;
  default:               func = 0; /* original PBE */
  }

  ss = X2S*x;

  f0 = kappa[func] + mu[func]*ss*ss;
  *f = 1.0 + kappa[func]*(1.0 - kappa[func]/f0);

  if(dfdx==NULL && d2fdx2==NULL) return; /* nothing else to do */

  df0 = 2.0*ss*mu[func];

  if(dfdx!=NULL){
    *dfdx  = X2S*kappa[func]*kappa[func]*df0/(f0*f0);
    *ldfdx = X2S*X2S*mu[func];
  }

  if(d2fdx2==NULL) return; /* nothing else to do */

  d2f0 = 2.0*mu[func];
  *d2fdx2 = X2S*X2S*kappa[func]*kappa[func]/(f0*f0)*(d2f0 - 2.0*df0*df0/f0);
}

void 
XC(gga_x_pbe_enhance)(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  func(p, x, f, dfdx, ldfdx, d2fdx2);
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pbe) = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbe_r) = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Revised PBE from Zhang & Yang",
  XC_FAMILY_GGA,
  "Y Zhang and W Yang, Phys. Rev. Lett 80, 890 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_pbe_sol) = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  "JP Perdew, et al, Phys. Rev. Lett. 100, 136406 (2008)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_x_xpbe) = {
  XC_GGA_X_XPBE,
  XC_EXCHANGE,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  "X Xu and WA Goddard III, J. Chem. Phys. 121, 4068 (2004)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL, NULL,
  work_gga_x
};
