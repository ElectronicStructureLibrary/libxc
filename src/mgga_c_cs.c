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

#define XC_MGGA_C_CS          72 /* Colle and Salvetti */

/*
    [1] Eq. (15) in http://dx.doi.org/10.1103/PhysRevB.37.785
    [2] CS2 in http://www.molpro.net/info/2012.1/doc/manual/node192.html

  there is a gamma(r) in [1] absent in [2]. This should be irrelevant
  for spin unpolarized. In any case, it seems that even in that case,
  libxc does not give the same as molpro, but I am unable to
  understand why...
*/

static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static FLOAT a = -0.04918, b = 0.132, c = 0.2533/RS_FACTOR, d = 0.349/RS_FACTOR;
  FLOAT cnst_283, opz, omz, opz2, omz2, opz13, omz13, opz23, omz23;
  FLOAT ta, tb, tw, aux, gamma, ff, num, den;
  FLOAT dgamma, dtwdz, dtwdxt, dtwdus0, dtwdus1, dauxdrs;
  FLOAT dnumdrs, dnumdz, dnumdxt, dnumdts0, dnumdts1, dnumdus0, dnumdus1, ddendrs;

  cnst_283 = 1.0/(4.0*M_CBRT2*M_CBRT2);

  opz = 1.0 + r->zeta;
  omz = 1.0 - r->zeta;

  opz2 = opz*opz;
  omz2 = omz*omz;

  opz13 = CBRT(opz); opz23 = opz13*opz13;
  omz13 = CBRT(omz); omz23 = omz13*omz13;

  ta    = r->ts[0] - r->us[0]/8.0;
  tb    = r->ts[1] - r->us[1]/8.0;
  tw    = r->xt*r->xt/8.0 - cnst_283*(r->us[0]*opz*opz23 + r->us[1]*omz*omz23)/4.0;

  aux   = EXP(-c*r->rs);
  gamma = 1.0 - r->zeta*r->zeta;

  ff    = cnst_283*(opz2*opz23*ta + omz2*omz23*tb) - tw;
  num   = (1.0 + 2.0*b*ff*aux)*gamma;
  den   = 1.0 + d*r->rs;

  r->f = a*num/den;

  if(r->order < 1) return;

  dtwdz   = -cnst_283*(5.0/3.0)*(r->us[0]*opz23 - r->us[1]*omz23)/4.0;
  dtwdxt  =  r->xt/4.0;
  dtwdus0 = -cnst_283*opz*opz23/4.0;
  dtwdus1 = -cnst_283*omz*omz23/4.0;

  dauxdrs = -c*aux;
  dgamma  = -2.0*r->zeta;

  dnumdrs  =  2.0*b*ff*dauxdrs*gamma;
  dnumdz   =  2.0*b*(cnst_283*(8.0/3.0)*(opz*opz23*ta - omz*omz23*tb) - dtwdz)*aux*gamma +
    (1.0 + 2.0*b*ff*aux)*dgamma;
  dnumdxt  = -2.0*b*dtwdxt*aux*gamma;
  dnumdts0 =  2.0*b*cnst_283*opz2*opz23*aux*gamma;
  dnumdts1 =  2.0*b*cnst_283*omz2*omz23*aux*gamma;
  dnumdus0 =  2.0*b*(-cnst_283*opz2*opz23/8.0 - dtwdus0)*aux*gamma;
  dnumdus1 =  2.0*b*(-cnst_283*omz2*omz23/8.0 - dtwdus1)*aux*gamma;

  ddendrs  = d;

  r->dfdrs    = a*(dnumdrs*den - num*ddendrs)/(den*den);
  r->dfdz     = a*dnumdz/den;
  r->dfdxt    = a*dnumdxt/den;
  r->dfdxs[0] = 0.0;
  r->dfdxs[1] = 0.0;
  r->dfdts[0] = a*dnumdts0/den;
  r->dfdts[1] = a*dnumdts1/den;
  r->dfdus[0] = a*dnumdus0/den;
  r->dfdus[1] = a*dnumdus1/den;

  if(r->order < 2) return;
}


#include "work_mgga_c.c"


const XC(func_info_type) XC(func_info_mgga_c_cs) = {
  XC_MGGA_C_CS,
  XC_CORRELATION,
  "Colle and Salvetti",
  XC_FAMILY_MGGA,
  {&xc_ref_Colle1975_329, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  NULL,
  NULL, NULL, NULL,
  work_mgga_c,
};

