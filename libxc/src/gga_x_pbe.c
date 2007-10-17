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

static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  static const double kappa[3] = {
    0.8040, /* original PBE */
    1.245,  /* PBE R */
    0.8040  /* PBE sol */
  };

  /* 
     the variable used in the original PBE paper is s = |grad n|/(2 k_f n)
     while here we use x = |grad n_s|/n_s^(4/3). Therefore, the value of mu
     we use here is he original mu multiplied by 

       1/2^(2/3) * 1/(4*(3 pi^2)^(2/3))

     where the first term comes from using the spin densities, and the second
     from the constants of k_f = (3 pi^2 n)^(1/3)
  */
  static const double mu[3] = {
    0.00361218645365094697,  /* PBE: mu = beta*pi^2/3, beta = 0.066725 (plus above mentioned constants) */
    0.00361218645365094697,  /* PBE rev: as PBE */
    0.00203151948716303243   /* PBE sol: 10/81 */
  };

  double dd;
  int func;

  switch(p->info->number){
  case XC_GGA_X_PBE_R:   func = 1; break;
  case XC_GGA_X_PBE_SOL: func = 2; break;
  default:               func = 0; /* original PBE */
  }

  dd     = 1.0/(kappa[func] + mu[func]*x*x);

  *f     = 1.0 + kappa[func]*(1.0 - kappa[func]*dd);
  *dfdx  = 2.0*x*mu[func]*kappa[func]*kappa[func]*dd*dd;
  *ldfdx = mu[func];
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_pbe = {
  XC_GGA_X_PBE,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

const xc_func_info_type func_info_gga_x_pbe_r = {
  XC_GGA_X_PBE_R,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)\n"
  "Y. Zhang and W. Yang, Phys. Rev. Lett 80, 890 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};

const xc_func_info_type func_info_gga_x_pbe_sol = {
  XC_GGA_X_PBE_SOL,
  XC_EXCHANGE,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "J.P.Perdew, K.Burke, and M.Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)\n"
  "J.P. Perdew, et al., arXiv:0707.2088v1",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
