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

#include "util.h"

#define XC_GGA_X_AK13  56 /* Armiento & Kuemmel 2013 */

static const double B1 =  1.74959015598863046792081721182; /* 3*muGE/5 + 8 pi/15 */
static const double B2 = -1.62613336586517367779736042170; /* muGE - B1 */

double XC(gga_ak13_get_asymptotic) (double homo)
{
  double Qx, aa, aa2, factor;

  Qx = SQRT(2.0)*B1/(3.0*CBRT(3.0*M_PI*M_PI));

  aa  = X_FACTOR_C*Qx;
  aa2 = aa*aa;

  factor = (homo < 0.0) ? -1.0 : 1.0;
    
  return (aa2/2.0)*(1.0 + factor*SQRT(1.0 - 4.0*homo/aa2));
}


#include "maple2c/gga_x_ak13.c"

#define func XC(gga_x_ak13_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_ak13) = {
  XC_GGA_X_AK13,
  XC_EXCHANGE,
  "Armiento & Kuemmel 2013",
  XC_FAMILY_GGA,
  {&xc_ref_Armiento2013_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  NULL, NULL, NULL,
  work_gga_x,
  NULL
};

