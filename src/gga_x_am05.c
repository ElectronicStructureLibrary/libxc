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

#include <config.h>
#if defined(HAVE_GSL)
#include <gsl/gsl_sf_lambert.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_AM05         120 /* Armiento & Mattsson 05 exchange                */

static void
lambert(FLOAT x, FLOAT *w, FLOAT *dw)
{
#if defined(HAVE_GSL)
  *w  = (FLOAT)gsl_sf_lambert_W0((double)x);
  if(x != 0.0)
    *dw = (*w)/(x*(1+(*w)));
  else
    *dw = 1.0;
#else
  fprintf(stderr, "Need the GSL library for this functional\n");
  exit(1);
#endif
}

static void /* Eq. (7) */
zeta_t(FLOAT s, FLOAT *f, FLOAT *df)
{
  FLOAT x, w, dw;

  x = POW(s, 1.5)/(2.0*sqrt(6.0));
  lambert(x, &w, &dw);

  *f  = POW(1.5*w, 2.0/3.0);
  if(s!=0)
    *df = POW(1.5, 1.0/6.0)*sqrt(s)*dw/(4.0*POW(w, 1.0/3.0));
  else
    *df = POW(6.0, 1.0/3.0)/4.0;
}

static void /* Eq. (3) */
zeta_tt(FLOAT s, FLOAT *f, FLOAT *df)
{
  FLOAT c;
  FLOAT z_t, z_t2, dz_t;
  zeta_t(s, &z_t, &dz_t);
  z_t2 = z_t*z_t;

  c   = POW(POW(4.0/3.0, 1.0/3.0) * 2.0*M_PI/3.0, 4);

  *f  = POW(c*z_t2 + z_t2*z_t2, 1.0/4.0);
  *df = POW(c*z_t2 + z_t2*z_t2, -3.0/4.0)*(2.0*c*z_t + 4.0*z_t2*z_t)*dz_t/4.0;
}

static void /* Eq. (10) */
n0_t(FLOAT s, FLOAT *f, FLOAT *df)
{
  FLOAT z_t, dz_t;

  zeta_t(s, &z_t, &dz_t);

  *f  = POW(z_t, 3.0/2.0)/(3.0*M_PI*M_PI*s*s*s);
  *df = (*f) * (3.0/2.0*dz_t/z_t - 3.0/s);
}

static void /* Eq. (3) */
ex_LDA(FLOAT n, FLOAT *f, FLOAT *df)
{
  FLOAT c;

  c = -3.0/(4.0*M_PI)*POW(3.0*M_PI*M_PI, 1.0/3.0);

  *f  = c*POW(n, 1.0/3.0);
  *df = c*POW(n, -2.0/3.0)/3.0;
}

static void /* Eq. (3) */
Fx_b(FLOAT s, FLOAT *f, FLOAT *df)
{
  FLOAT n0, dn0, z_tt, dz_tt, ex, dex;

  n0_t(s, &n0, &dn0);
  zeta_tt(s, &z_tt, &dz_tt);
  ex_LDA(n0, &ex, &dex);

  *f  = -1.0/(ex*2.0*z_tt);
  *df = 2.0*(*f)*(*f) * (dex*dn0*z_tt + ex*dz_tt);
}

static void /* Eq. (12) */
XX(FLOAT s, FLOAT *f, FLOAT *df)
{
  const FLOAT alpha = 2.804;
  FLOAT s2, f1;

  s2 = s*s;
  f1 = 1.0 + alpha*s2;

  *f  = 1.0 - alpha*s2/f1;
  *df = -2.0*alpha*s/(f1*f1);
}

static void /* Eq. (6) */
F_LAA(FLOAT s, FLOAT *f, FLOAT *df)
{
  const FLOAT c=0.7168;
  FLOAT s2, fxb, dfxb, f1, f2, df1, df2;

  s2 = s*s;
  Fx_b(s, &fxb, &dfxb);

  f1  = c*s2 + 1.0;
  df1 = 2.0*c*s;

  f2  = c*s2/fxb + 1.0;
  df2 = c*(2.0*s*fxb - dfxb*s2)/(fxb*fxb);

  *f  = f1/f2;
  *df = (df1*f2 - f1*df2)/(f2*f2);
}

static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  const FLOAT x2s     = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  double ss, xx, dxx, flaa, dflaa;
   
  if(x >= MIN_GRAD){
    ss  = x2s*x;

    XX(ss, &xx, &dxx);
    F_LAA(ss, &flaa, &dflaa);

    *f     = xx + (1.0 - xx)*flaa;
    *dfdx  = dxx*(1.0 - flaa) + dflaa*(1.0 - xx);
  }else{
    *f    = 1.0;
  }

  *ldfdx = -2.804; /* -alpha?? */

  *dfdx  *= x2s;
  *ldfdx *= x2s*x2s;
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_am05) = {
  XC_GGA_X_AM05,
  XC_EXCHANGE,
  "Armiento & Mattsson 05",
  XC_FAMILY_GGA,
  "R Armiento and AE Mattsson, Phys. Rev. B 72, 085108 (2005)\n"
  "AE Mattsson, R Armiento, J Paier, G Kresse, JM Wills, and TR Mattsson, J. Chem. Phys. 128, 084714 (2008).",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
