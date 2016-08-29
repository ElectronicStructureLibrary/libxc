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

#define XC_MGGA_XC_ZLP          42 /* Zhao, Levy & Parr, Eq. (21) */

static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static FLOAT cc = 0.828432*RS_FACTOR, dd = 2.15509e-2*RS_FACTOR, kk = 2.047107e-3*RS_FACTOR;
  FLOAT cnst_253, tw, aux1, aux2, aux3, opz, omz, opz2, omz2, opz13, omz13, opz23, omz23;
  FLOAT dtwdz, dtwdxt, dtwdus0, dtwdus1, daux2, daux3;

  cnst_253 = 1.0/(2.0*M_CBRT2*M_CBRT2);

  opz = 1.0 + r->zeta;
  omz = 1.0 - r->zeta;

  opz2 = opz*opz;
  omz2 = omz*omz;

  opz13 = CBRT(opz); opz23 = opz13*opz13;
  omz13 = CBRT(omz); omz23 = omz13*omz13;

  tw    = (r->xt*r->xt - cnst_253*(r->us[0]*opz*opz23 + r->us[1]*omz*omz23))/8.0; 

  aux1 = -(cc + dd*tw);
  aux2 = LOG(1.0 + r->rs/kk);
  aux3 = (1.0 - kk*aux2/r->rs)/r->rs;

  r->f = aux1*aux3;

  if(r->order < 1) return;

  dtwdz   = -cnst_253*(5.0/3.0)*(r->us[0]*opz23 - r->us[1]*omz23)/8.0;
  dtwdxt  =  r->xt/4.0;
  dtwdus0 = -cnst_253*opz*opz23/8.0;
  dtwdus1 = -cnst_253*omz*omz23/8.0;

  daux2 = 1.0/(r->rs + kk);
  daux3 = -(r->rs - 2.0*kk*aux2 + kk*r->rs*daux2)/(r->rs*r->rs*r->rs);

  r->dfdrs    = aux1*daux3;
  r->dfdz     = -dd*dtwdz*aux3;
  r->dfdxt    = -dd*dtwdxt*aux3;
  r->dfdxs[0] = 0.0;
  r->dfdxs[1] = 0.0;
  r->dfdts[0] = 0.0;
  r->dfdts[1] = 0.0;
  r->dfdus[0] = -dd*dtwdus0*aux3;
  r->dfdus[1] = -dd*dtwdus1*aux3;

  if(r->order < 2) return;
}


#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_xc_zlp) = {
  XC_MGGA_XC_ZLP,
  XC_EXCHANGE_CORRELATION,
  "Zhao, Levy & Parr, Eq. (21)",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao1993_918, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  NULL,
  NULL, NULL, NULL,
  work_mgga_c,
};

