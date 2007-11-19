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

#define XC_GGA_X_PW91         109 /* Perdew & Wang 91 */

static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  const double x2s = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  const double aa = 0.19645, bb = 7.7956, cc = 0.2743, dd=-0.1508, ff=0.004, alpha=100.0;
  double ss, ss2, ss4;
  double f1, f2, f3, f4;

  ss  = x2s*x;
  ss2 = ss*ss;
  ss4 = ss2*ss2;

  f1 = dd*exp(-alpha*ss2);
  f2 = aa*asinh(bb*ss);
  f3 = (cc + f1)*ss2 - ff*ss4;
  f4 = 1.0 + ss*f2 + ff*ss4;

  *f     = 1.0 + f3/f4;
  *dfdx  = (2.0*ss*(cc + f1*(1.0 - alpha*ss2) - 2.0*ff*ss2)*f4 - 
	    f3*(f2 + ss*(aa*bb/sqrt(1.0 + bb*bb*ss2) + 4.0*ff*ss2)))/(f4*f4);
  *dfdx *= x2s;
  *ldfdx = x2s*x2s*(cc + dd);
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_pw91 = {
  XC_GGA_X_PW91,
  XC_EXCHANGE,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, and C Fiolhais, Phys. Rev. B 46, 6671 (1992)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
