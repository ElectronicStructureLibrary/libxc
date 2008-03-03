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

#define XC_GGA_X_FT97_A       114 /* Filatov & Thiel 97 (version A) */
#define XC_GGA_X_FT97_B       115 /* Filatov & Thiel 97 (version B) */

static inline void
func(xc_gga_type *p, FLOAT x, FLOAT sigma, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *vsigma)
{
  static const FLOAT 
    beta0 = 0.002913644, beta1 = 0.0009474169, beta2 = 6255746.320201; /* beta2 = 2501.149^2 ?? (Eq. (16a) */

  FLOAT x2, beta, dbetadsigma, f1, f2, f3;
  int func;

  switch(p->info->number){
  case XC_GGA_X_FT97_B: func = 1; break;
  default:              func = 0; /*  XC_GGA_X_FT97_A */
  }

  if(func == 0){
    beta = 0.00293;
    dbetadsigma = 0.0;
  }else{
    f1   = beta2 + sigma;
    beta = beta0 + beta1*sigma/f1;
    dbetadsigma = beta1*beta2/(f1*f1);
  }

  x2 = x*x;
  f2 = beta*asinh(x2);
  f3 = sqrt(1.0 + 9.0*x2*f2*f2);
  *f = 1.0 + beta/X_FACTOR_C*x2/f3;
 
  *dfdx = beta/X_FACTOR_C*2.0*x*( f3 - 9.0/2.0*x2/f3*(f2*f2 + 2.0*x2*beta*f2/sqrt(1 + x2*x2)) )/(f3*f3);
  *ldfdx= beta0/X_FACTOR_C;

  *vsigma = dbetadsigma*x2/(f3*X_FACTOR_C)*(1.0 - 9.0*x2*f2*f2/(f3*f3));
}

#define HEADER 2
#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_ft97_a = {
  XC_GGA_X_FT97_A,
  XC_EXCHANGE,
  "Filatov & Thiel 97 (version A)",
  XC_FAMILY_GGA,
  "M Filatov and W Thiel, Mol. Phys 91, 847 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL, 
  work_gga_x
};

const xc_func_info_type func_info_gga_x_ft97_b = {
  XC_GGA_X_FT97_B,
  XC_EXCHANGE,
  "Filatov & Thiel 97 (version B)",
  XC_FAMILY_GGA,
  "M Filatov and W Thiel, Mol. Phys 91, 847 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL, 
  work_gga_x
};
