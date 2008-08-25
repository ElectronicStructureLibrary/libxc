/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

#include <stdlib.h>
#include <assert.h>
#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_LTA          201 /* Local tau approximation of Ernzerhof & Scuseria */

static void 
func(const XC(mgga_type) *p, FLOAT x, FLOAT tau, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *dfdtau, 
     FLOAT *d2fdx2, FLOAT *d2fdxtau, FLOAT *d2fdtau2)
{
  const FLOAT a1 = 0.430075922439080216009; /* POW(10.0/(3.0*POW(3.0*M_PI*M_PI, 2.0/3.0)), 4.0/5.0) */

  *f = a1*POW(tau, 4.0/5.0);
  if(dfdx==NULL && d2fdx2==NULL) return; /* nothing else to do */

  *dfdx   = 0.0;
  *ldfdx  = 0.0;
  *dfdtau = a1*4.0/5.0*POW(tau, -1.0/5.0);

  if(d2fdx2==NULL) return; /* nothing else to do */

  *d2fdx2   = 0.0;
  *d2fdxtau = 0.0;
  *d2fdtau2 = -a1*4.0/25.0*POW(tau, -6.0/5.0);
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_lta) = {
  XC_MGGA_X_LTA,
  XC_EXCHANGE,
  "Local tau approximation",
  XC_FAMILY_MGGA,
  "M Ernzerhof and G Scuseria, J. Chem. Phys. 111, 911 (1999)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
