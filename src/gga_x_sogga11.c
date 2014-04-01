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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_SOGGA11        151 /* Second-order generalized gradient approximation 2011 */
#define XC_HYB_GGA_X_SOGGA11_X  426 /* Hybrid based on SOGGA11 form */

static void 
gga_x_sogga11_init(XC(func_type) *p)
{
  switch(p->info->number){
  case XC_GGA_X_SOGGA11:
    p->func = 0;
    break;
  case XC_HYB_GGA_X_SOGGA11_X:
    p->func = 1;
    p->cam_alpha = 0.4015;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_sogga11\n");
    exit(1);
  }
}

void XC(gga_x_sogga11_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  const FLOAT kappa = 0.552;
  const FLOAT mu = 10.0/81.0;
  const FLOAT alpha = mu*X2S*X2S/kappa;
  const FLOAT aa[][6] = {
    {0.50000, -2.95535,  15.7974, -91.1804,  96.2030, 0.18683},   /* SOGGA11   */
    {0.50000,  5.37406, -5.94160,  12.7962, -18.8521, 8.78551}    /* SOGGA11-X */
  };
  const FLOAT bb[][6] = {
    {0.50000,  3.50743, -12.9523,  49.7870, -33.2545, -11.1396},  /* SOGGA11   */
    {0.50000, -4.82197,   5.40713, -4.10014, -6.27393,  6.62678}  /* SOGGA11-X */
  };
    
  FLOAT f0, df0, d2f0, d3f0, den0, den1, t0, dt0, d2t0, d3t0, t1, dt1, d2t1, d3t1, f1, df1, d2f1, d3f1;

  den0 = -1.0/(1.0 + alpha*x*x);
  f0   =  1.0 + den0;
  den1 = -EXP(-alpha*x*x);
  f1   =  1.0 + den1;

  t0 = aa[p->func][0] + f0*(aa[p->func][1] + f0*(aa[p->func][2] + f0*(aa[p->func][3] + f0*(aa[p->func][4] + f0*aa[p->func][5]))));
  t1 = bb[p->func][0] + f1*(bb[p->func][1] + f1*(bb[p->func][2] + f1*(bb[p->func][3] + f1*(bb[p->func][4] + f1*bb[p->func][5]))));

  *f  = t0;
  *f += t1;

  if(order < 1) return;

  df0 =  2.0*alpha*x*den0*den0;
  df1 = -2.0*alpha*x*den1;

  dt0  = aa[p->func][1] + f0*(2.0*aa[p->func][2] + f0*(3.0*aa[p->func][3] + f0*(4.0*aa[p->func][4] + f0*5.0*aa[p->func][5])));
  dt1  = bb[p->func][1] + f1*(2.0*bb[p->func][2] + f1*(3.0*bb[p->func][3] + f1*(4.0*bb[p->func][4] + f1*5.0*bb[p->func][5])));

  *dfdx  = dt0*df0;
  *dfdx += dt1*df1;

  if(order < 2) return;

  d2f0 = 2.0*alpha*(3.0*alpha*x*x - 1.0)*den0*den0*den0;
  d2f1 = 2.0*alpha*(2.0*alpha*x*x - 1.0)*den1;

  d2t0 = 2.0*aa[p->func][2] + f0*(3.0*2.0*aa[p->func][3] + f0*(4.0*3.0*aa[p->func][4] + f0*5.0*4.0*aa[p->func][5]));
  d2t1 = 2.0*bb[p->func][2] + f1*(3.0*2.0*bb[p->func][3] + f1*(4.0*3.0*bb[p->func][4] + f1*5.0*4.0*bb[p->func][5]));

  *d2fdx2  = dt0*d2f0 + df0*df0*d2t0;
  *d2fdx2 += dt1*d2f1 + df1*df1*d2t1;

  if(order < 3) return;

  d3f0 = 24.0*alpha*alpha*x*(alpha*x*x - 1.0)*den0*den0*den0*den0;
  d3f1 = -4.0*alpha*alpha*x*(2.0*alpha*x*x - 3.0)*den1;

  d3t0 = 3.0*2.0*aa[p->func][3] + f0*(4.0*3.0*2.0*aa[p->func][4] + f0*5.0*4.0*3.0*aa[p->func][5]);
  d3t1 = 3.0*2.0*bb[p->func][3] + f1*(4.0*3.0*2.0*bb[p->func][4] + f1*5.0*4.0*3.0*bb[p->func][5]);

  *d3fdx3  = 3.0*df0*d2f0*d2t0 + dt0*d3f0 + df0*df0*df0*d3t0;
  *d3fdx3 += 3.0*df1*d2f1*d2t1 + dt1*d3f1 + df1*df1*df1*d3t1;
}


#define func XC(gga_x_sogga11_enhance)
#include "work_gga_x.c"


const XC(func_info_type) XC(func_info_gga_x_sogga11) = {
  XC_GGA_X_SOGGA11,
  XC_EXCHANGE,
  "Second-order generalized gradient approximation 2011",
  XC_FAMILY_GGA,
  "R Peverati, Y Zhao, and DG Truhlar, J. Phys. Chem. Lett. 2, 1991-1997 (2011)\n"
  "http://comp.chem.umn.edu/mfm/index.html",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-31, 1e-32, 0.0, 1e-32,
  gga_x_sogga11_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_x_sogga11_x) = {
  XC_HYB_GGA_X_SOGGA11_X,
  XC_EXCHANGE,
  "Hybrid based on SOGGA11 form",
  XC_FAMILY_HYB_GGA,
  "R Peverati and DG Truhlar, J. Chem. Phys. 135, 191102 (2011)\n"
  "http://comp.chem.umn.edu/mfm/index.html",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-31, 1e-32, 0.0, 1e-32,
  gga_x_sogga11_init, 
  NULL, NULL,
  work_gga_x,
  NULL
};
