/*
 Copyright (C) 2006-2009 M.A.L. Marques

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "util.h"

#define XC_MGGA_X_BR89         206 /* Becke-Roussel 89  */

static FLOAT br89_gamma = 0.8;

/* This code follows the inversion done in the PINY_MD package */

FLOAT inline 
br_newt_raph(FLOAT a, FLOAT tol, int *ierr)
{
  int count;
  double x, f;
  static int max_iter = 100;

   *ierr = 1;
   if(a == 0.0)
     return 0.0;
   
   /* starting point */
   if(a < 0.0)
     x = 2.0 - tol;
   else{ /* a > 0 */
     x = -0.375*(log(a) - 5.0);
     if(x <= 2.0) x = 2.0 + tol;
   }

   count = 0;
   do {
     double arg, eee, xm2, fp;

     xm2 = x - 2.0;
     arg = 2.0*x/3.0;
     eee = exp(-arg);

     f  = x*eee - a*xm2;
     fp = eee*(1.0 - 2.0/3.0*x) - a;

     x -= f/fp;
     x  = fabs(x);

     count ++;
   }while((fabs(f) > tol) && (count < max_iter));

   if(count == max_iter) *ierr=0; 
   return x;
}


FLOAT inline
br_bisect(FLOAT a, FLOAT tol, int *ierr)
{
   int count;
   FLOAT f, x, x1, x2;
   static int max_iter = 500;

   *ierr = 1;
   if(a == 0.0)
     return 0.0;

   /* starting interval */
   if(a > 0.0) {
     x1 = 2.0 + tol;
     x2 = 10.0;
   }else{
     x2 = 2.0 - tol;
     x1 = 0.0;
   }

   /* bisection */
   count = 0;
   do{
     FLOAT arg, eee, xm2;
     x   = 0.5*(x1 + x2);
     xm2 = x - 2.0;
     arg = 2.0*x/3.0;
     eee = exp(-arg);
     f   = x*eee - a*xm2;

     if(f > 0.0) x1 = x;
     if(f < 0.0) x2 = x;

     count++;
   }while((fabs(f) > tol)  && (count < max_iter));

   if(count == max_iter) *ierr=0; 
   return x;
}

FLOAT XC(mgga_x_br89_get_x)(FLOAT Q)
{
  FLOAT rhs, br_x, tol;
  int ierr;

#if SINGLE_PRECISION
  tol = 1e-6;
#else
  tol = 5e-12;
#endif

  /* build right-hand side of the non-linear equation 
     Remember we use a different definition of tau */
  rhs = 2.0/3.0*POW(M_PI, 2.0/3.0)/Q;

  br_x = br_newt_raph(rhs, tol, &ierr);
  if(ierr == 0 || isnan(br_x) != 0)
    br_x = br_bisect(rhs, tol, &ierr);
  if(ierr == 0){
    fprintf(stderr, "Warning: Convergence not reached in Becke-Roussel functional\n");
  }

  return br_x;
}

static void 
func(const XC(mgga_type) *pt, FLOAT x, FLOAT t, FLOAT u, int order,
     FLOAT *f, FLOAT *dfdx, FLOAT *dfdt, FLOAT *dfdu,
     FLOAT *d2fdx2, FLOAT *d2fdxt, FLOAT *d2fdt2)
{
  FLOAT Q, br_x, dfdbx, dxdu, ff, dff;
  FLOAT cnst, exp1, exp2;
 
  Q  = (u - 2.0*br89_gamma*t + 0.5*br89_gamma*x*x)/6.0;
  br_x = XC(mgga_x_br89_get_x)(Q);

  /* we have also to include the factor 1/2 from Eq. (9) */
  cnst = POW(M_PI, 1.0/3.0)/X_FACTOR_C;
  exp1 = exp(br_x/3.0);
  exp2 = exp(-br_x);

  *f = cnst*exp1*(1.0 - exp2*(1.0 + br_x/2.0))/br_x;

  if(order < 1) return;

  dfdbx = cnst*(3.0 + br_x*(br_x + 2.0) + (br_x - 3.0)/exp2) /
    (3.0*exp1*exp1*br_x*br_x);
  
  ff  = br_x*exp(-2.0/3.0*br_x)/(br_x - 2);
  dff = -2.0/3.0 + 1.0/br_x - 1.0/(br_x - 2.0); /* dff / ff */

  dxdu  = -1.0/(6.0*Q*dff);

  *dfdx =    x*br89_gamma*dfdbx*dxdu;
  *dfdt = -2.0*br89_gamma*dfdbx*dxdu;
  *dfdu =                 dfdbx*dxdu;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_br89) = {
  XC_MGGA_X_BR89,
  XC_EXCHANGE,
  "Becke-Roussel 89",
  XC_FAMILY_MGGA,
  "AD Becke and MR Roussel, Phys. Rev. A 39, 3761 (1989)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
