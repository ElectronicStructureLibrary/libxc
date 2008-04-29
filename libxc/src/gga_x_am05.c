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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_AM05         120 /* Armiento & Mattsson 05 exchange                */

static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  const FLOAT am05_c      = 0.7168;
  const FLOAT am05_alpha  = 2.804;

  const FLOAT z_tt_factor = POW(POW(4.0/3.0, 1.0/3.0) * 2.0*M_PI/3.0, 4);

  FLOAT ss, ss2, lam_x, ww;
  FLOAT z_t, z_t2, z_tt, fx_b, xx, flaa_1, flaa_2, flaa;

  if(x < MIN_GRAD){
    *f    = 1.0;
    if(dfdx != NULL) *ldfdx = -am05_alpha*X2S*X2S;
    return;
  }

  ss  = X2S*x;
  ss2 = ss*ss;

  lam_x  = POW(ss, 1.5)/(2.0*sqrt(6.0));
  ww     = (FLOAT)lambert_w((double)lam_x);

  z_t    = POW(1.5*ww, 2.0/3.0);
  z_t2   = z_t*z_t;

  /* This is equal to sqrt(t_zeta) * tt_zeta} of the JCP*/
  z_tt   = z_t * POW(z_tt_factor + z_t2, 1.0/4.0);

  /* note that there is a factor of 2 missing in the JCP */
  fx_b   = M_PI/3.0*ss/z_tt;

  xx     = 1.0/(1.0 + am05_alpha*ss2);

  flaa_1 = am05_c*ss2 + 1.0;
  flaa_2 = am05_c*ss2/fx_b + 1.0;
  flaa   = flaa_1/flaa_2;

  *f     = xx + (1.0 - xx)*flaa;

  if(dfdx != NULL){
    FLOAT dww, dz_t, dz_tt, dfx_b, dxx, dflaa_1, dflaa_2, dflaa;

    dww     = ww/(lam_x*(1.0 + ww));
    dz_t    = POW(1.5, 1.0/6.0)*sqrt(ss)*dww/(4.0*POW(ww, 1.0/3.0));
    dz_tt   = POW(z_tt_factor + z_t2, -3.0/4.0)*(z_tt_factor + 3.0*z_t2/2.0)*dz_t;
    dfx_b   = M_PI/3.0*(z_tt - ss*dz_tt)/(z_tt*z_tt);

    dxx     = -2.0*am05_alpha*ss * xx*xx;
    dflaa_1 = 2.0*am05_c*ss;
    dflaa_2 = am05_c*(2.0*ss*fx_b - dfx_b*ss2)/(fx_b*fx_b);
    dflaa   = (dflaa_1*flaa_2 - flaa_1*dflaa_2)/(flaa_2*flaa_2);

    *dfdx  = dxx*(1.0 - flaa) + dflaa*(1.0 - xx);
    *ldfdx = -am05_alpha; /* -alpha?? */

    *dfdx  *= X2S;
    *ldfdx *= X2S*X2S;
  }

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
