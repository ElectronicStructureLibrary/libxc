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

#define XC_MGGA_C_PKZB          239 /* Perdew, Kurth, Zupan, and Blaha */

static void 
mgga_c_pkzb_init(XC(func_type) *p)
{
  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_C_PBE, p->nspin);
}


static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static FLOAT C = 0.53;
  static const FLOAT tmin = 0.5e-10;

  XC(gga_work_c_t) PBE[3];
  FLOAT opz, omz, opz13, omz13, opz23, omz23, taut, xtot, dd, dd2;
  FLOAT dtautdz, dtautts[2], dxtotdz, dxtotdxs[2];
  FLOAT ddddz, ddddxs[2], ddddts[2];
  int is;

  /* first we get the parallel and perpendicular PBE */
  XC(pbe_c_stoll) (pt->func_aux[0], r, PBE);

  opz = 1.0 + r->zeta;
  omz = 1.0 - r->zeta;

  opz13 = CBRT(opz);
  omz13 = CBRT(omz);

  opz23 = opz13*opz13;
  omz23 = omz13*omz13;

  /* initialize to zero */
  r->f = 0.0;
  if(r->order >= 1){
    r->dfdrs = r->dfdz = r->dfdxs[0] = r->dfdxs[1] = r->dfdxt = 0.0;
    r->dfdus[0] = r->dfdus[1] = r->dfdts[0] = r->dfdts[1] = 0.0;
  }

  for(is = 0; is < 2; is++){
    dd  = (r->ts[is] > tmin) ? r->xs[is]*r->xs[is]/(8.0*r->ts[is]) : 0.0;
    dd2 = dd*dd;

    r->f += -(1.0 + C)*dd*dd*PBE[is].f;

    if(r->order < 1) continue;

    if(r->ts[is] > tmin){
      ddddxs[is] = 2.0*r->xs[is]/(8.0*r->ts[is]);
      ddddts[is] = -r->xs[is]*r->xs[is]/(8.0*r->ts[is]*r->ts[is]);
    }else
      ddddxs[is] = ddddts[is] = 0.0;

    r->dfdrs     += -(1.0 + C)*dd2*PBE[is].dfdrs;
    r->dfdz      += -(1.0 + C)*dd2*PBE[is].dfdz;
    r->dfdxs[is] += -(1.0 + C)*dd*(2.0*ddddxs[is]*PBE[is].f + dd*PBE[is].dfdxs[is]);
    r->dfdts[is] += -(1.0 + C)*2.0*dd*ddddts[is]*PBE[is].f;
  }

  taut = r->ts[0]*opz*opz23 + r->ts[1]*omz*omz23;
  xtot = r->xs[0]*r->xs[0]*opz*opz23 + r->xs[1]*r->xs[1]*omz*omz23;

  dd = (taut > tmin) ? xtot/(8.0*taut) : 0.0;
  dd2 = dd*dd;
  
  r->f += (1.0 + C*dd2)*PBE[2].f;

  if(r->order < 1) return;

  if(taut > tmin){
    dtautdz     = 5.0/3.0 * (r->ts[0]*opz23 - r->ts[1]*omz23);
    dtautts[0]  = opz*opz23;
    dtautts[1]  = omz*omz23;

    dxtotdz     = 5.0/3.0 * (r->xs[0]*r->xs[0]*opz23 - r->xs[1]*r->xs[1]*omz23);
    dxtotdxs[0] = 2.0*r->xs[0]*opz*opz23;
    dxtotdxs[1] = 2.0*r->xs[1]*omz*omz23;

    ddddz     = (dxtotdz*taut - xtot*dtautdz)/(8.0*taut*taut);
    ddddxs[0] = dxtotdxs[0]/(8.0*taut);
    ddddxs[1] = dxtotdxs[1]/(8.0*taut);
    ddddts[0] = -xtot*dtautts[0]/(8.0*taut*taut);
    ddddts[1] = -xtot*dtautts[1]/(8.0*taut*taut);
  }else{
    ddddz = ddddxs[0] = ddddxs[1] = ddddts[0] = ddddts[1] = 0.0;
  }

  r->dfdrs    += (1.0 + C*dd2)*PBE[2].dfdrs;
  r->dfdz     += (1.0 + C*dd2)*PBE[2].dfdz + 2.0*C*dd*ddddz*PBE[2].f;
  r->dfdxt     = (1.0 + C*dd2)*PBE[2].dfdxt;
  r->dfdxs[0] += (1.0 + C*dd2)*PBE[2].dfdxs[0] + 2.0*C*dd*ddddxs[0]*PBE[2].f;
  r->dfdxs[1] += (1.0 + C*dd2)*PBE[2].dfdxs[1] + 2.0*C*dd*ddddxs[1]*PBE[2].f;
  r->dfdts[0] += 2.0*C*dd*ddddts[0]*PBE[2].f;
  r->dfdts[1] += 2.0*C*dd*ddddts[1]*PBE[2].f;

  if(r->order < 2) return;

}


#include "work_mgga_c.c"


XC(func_info_type) XC(func_info_mgga_c_pkzb) = {
  XC_MGGA_C_PKZB,
  XC_CORRELATION,
  "Perdew, Kurth, Zupan, and Blaha",
  XC_FAMILY_MGGA,
  "JP Perdew, S Kurth, A Zupan, and P. Blaha, Phys. Rev. Lett. 82, 2544-2547 (1999)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_c_pkzb_init,
  NULL, NULL, NULL,
  work_mgga_c,
};
