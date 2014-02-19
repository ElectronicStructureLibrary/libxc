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
#include <string.h>
#include <assert.h>

#include "util.h"

#define XC_GGA_C_Q2D          47 /* Chiodo et al  */

static void gga_c_q2d_init(XC(func_type) *p)
{
  p->n_func_aux  = 2;
  p->func_aux    = (XC(func_type) **) malloc(2*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));
  p->func_aux[1] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_2D_AMGB, p->nspin);
  XC(func_init)(p->func_aux[1], XC_GGA_C_PBE, p->nspin);
}


inline void 
XC(gga_c_q2d_func) (const XC(func_type) *p, XC(gga_work_c_t) *r)
{
  const FLOAT rs2D_factor = 1.704, dd = 1e6;

  FLOAT rs2D, phi, tconv, auxp, auxm, t, t2, t4, t6, num, den, fac;
  FLOAT drs2Ddrs, drs2Ddxt, d2rs2Ddrsxt, d2rs2Ddxt2;
  FLOAT dphidz, dtdrs, dtdxt, dtdphi, dtdz, dnumdt, ddendt, dfacdt;
  FLOAT d2phidz2, d2tdrs2, d2tdrsxt, d2tdphi2, d2tdrsphi, d2tdxtphi, d2tdz2, d2tdrsz, d2tdzxt, d2numdt2, d2dendt2, d2facdt2;

  XC(lda_work_t) ldaw;
  XC(gga_work_c_t) ggaw;

  /* we start by getting the 2D LDA */
  rs2D = rs2D_factor*r->rs*SQRT(X2S*r->xt)/RS_FACTOR;

  ldaw.order = r->order;
  
  ldaw.rs[0] = SQRT(rs2D);
  ldaw.rs[1] = rs2D;
  ldaw.rs[2] = rs2D*rs2D;
  ldaw.zeta  = r->zeta;

  XC(lda_c_2d_amgb_func)(p->func_aux[0], &ldaw);

  /* now we get the PBE */
  memcpy(&ggaw, r, sizeof(XC(gga_work_c_t)));
  
  XC(gga_c_pbe_func)(p->func_aux[1], &ggaw);

  /* now comes the interpolation between them */
  tconv = 4.0*M_CBRT2;

  auxp = CBRT(1.0 + r->zeta);
  auxm = CBRT(1.0 - r->zeta);

  phi  = 0.5*(auxp*auxp + auxm*auxm);
  t    = r->xt/(tconv*phi*SQRT(r->rs));

  t2 = t*t;
  t4 = t2*t2;
  t6 = t4*t2;

  num = t4*(1.0 + t2);
  den = dd + t6;
  fac = num/den;

  r->f = ggaw.f + fac*(-ggaw.f + ldaw.zk);

  if(r->order < 1) return;

  drs2Ddrs = rs2D/r->rs;
  drs2Ddxt = rs2D/(2.0*r->xt);

  dphidz = 0.0;
  if(auxp > p->info->min_zeta) dphidz += 1/auxp;
  if(auxm > p->info->min_zeta) dphidz -= 1/auxm;
  dphidz *= 1.0/3.0;

  dtdrs  = -t/(2.0*r->rs);
  dtdxt  =  t/r->xt;
  dtdphi = -t/phi;
  dtdz   = dtdphi*dphidz;

  dnumdt = t*t2*(4.0 + 6.0*t2);
  ddendt = 6.0*t*t4;
  dfacdt = DFRACTION(num, dnumdt, den, ddendt);

  r->dfdrs    = ggaw.dfdrs + fac*(-ggaw.dfdrs + ldaw.dedrs*drs2Ddrs) + dfacdt*dtdrs*(-ggaw.f + ldaw.zk);
  r->dfdz     = ggaw.dfdz  + fac*(-ggaw.dfdz  + ldaw.dedz) + dfacdt*dtdz*(-ggaw.f + ldaw.zk);
  r->dfdxt    = ggaw.dfdxt + fac*(-ggaw.dfdxt + ldaw.dedrs*drs2Ddxt) + dfacdt*dtdxt*(-ggaw.f + ldaw.zk);
  r->dfdxs[0] = (1.0 - fac)*ggaw.dfdxs[0];
  r->dfdxs[1] = (1.0 - fac)*ggaw.dfdxs[1];


  if(r->order < 2) return;

  d2rs2Ddrsxt = drs2Ddxt/r->rs;
  d2rs2Ddxt2  = -drs2Ddxt/(2.0*r->xt);

  d2phidz2 = 0.0;
  if(auxp > p->info->min_zeta) d2phidz2 += 1.0/((1.0 + r->zeta)*auxp);
  if(auxm > p->info->min_zeta) d2phidz2 += 1.0/((1.0 - r->zeta)*auxm);
  d2phidz2 *= -1.0/9.0;

  d2tdrs2   =  -3.0*dtdrs/(2.0*r->rs);
  d2tdrsxt  =  dtdrs/r->xt;

  d2tdphi2  = -2.0*dtdphi/phi;
  d2tdrsphi = -dtdrs/phi;
  d2tdxtphi =  dtdphi/r->xt;

  d2tdz2    =  d2tdphi2*dphidz*dphidz + dtdphi*d2phidz2;
  d2tdrsz   =  d2tdrsphi*dphidz;
  d2tdzxt   =  d2tdxtphi*dphidz;

  d2numdt2 = t2*(12.0 + 30.0*t2);
  d2dendt2 = 30.0*t4;
  d2facdt2 = D2FRACTION(num, dnumdt, d2numdt2, den, ddendt, d2dendt2);

  r->d2fdrs2     = ggaw.d2fdrs2 + fac*(-ggaw.d2fdrs2 + ldaw.d2edrs2*drs2Ddrs*drs2Ddrs) + 
    2.0*dfacdt*dtdrs*(-ggaw.dfdrs + ldaw.dedrs*drs2Ddrs) + (d2facdt2*dtdrs*dtdrs + dfacdt*d2tdrs2)*(-ggaw.f + ldaw.zk);
  
  r->d2fdrsz     = ggaw.d2fdrsz + fac*(-ggaw.d2fdrsz + ldaw.d2edrsz*drs2Ddrs) + 
    dfacdt*dtdrs*(-ggaw.dfdz + ldaw.dedz) + dfacdt*dtdz*(-ggaw.dfdrs + ldaw.dedrs*drs2Ddrs) +
    (d2facdt2*dtdrs*dtdz + dfacdt*d2tdrsz)*(-ggaw.f + ldaw.zk);

  r->d2fdrsxt    =  ggaw.d2fdrsxt + fac*(-ggaw.d2fdrsxt + ldaw.d2edrs2*drs2Ddrs*drs2Ddxt + ldaw.dedrs*d2rs2Ddrsxt) + 
    dfacdt*dtdrs*(-ggaw.dfdxt + ldaw.dedrs*drs2Ddxt) + dfacdt*dtdxt*(-ggaw.dfdrs + ldaw.dedrs*drs2Ddrs) +
    (d2facdt2*dtdrs*dtdxt + dfacdt*d2tdrsxt)*(-ggaw.f + ldaw.zk);

  r->d2fdrsxs[0] = 0.0;
  r->d2fdrsxs[1] = 0.0;

  r->d2fdz2      = ggaw.d2fdz2 + fac*(-ggaw.d2fdz2 + ldaw.d2edz2) + 
    2.0*dfacdt*dtdz*(-ggaw.dfdz + ldaw.dedz) + (d2facdt2*dtdz*dtdz + dfacdt*d2tdz2)*(-ggaw.f + ldaw.zk);

  r->d2fdzxt     = ggaw.d2fdzxt + fac*(-ggaw.d2fdzxt + ldaw.d2edrsz*drs2Ddxt) + 
    dfacdt*dtdxt*(-ggaw.dfdz + ldaw.dedz) + dfacdt*dtdz*(-ggaw.dfdxt + ldaw.dedrs*drs2Ddxt) +
    (d2facdt2*dtdxt*dtdz + dfacdt*d2tdzxt)*(-ggaw.f + ldaw.zk);

  r->d2fdzxs[0]  = 0.0;
  r->d2fdzxs[1]  = 0.0;

  r->d2fdxt2     = ggaw.d2fdxt2 + fac*(-ggaw.d2fdxt2 + ldaw.d2edrs2*drs2Ddxt*drs2Ddxt + ldaw.dedrs*d2rs2Ddxt2) + 
    2.0*dfacdt*dtdxt*(-ggaw.dfdxt + ldaw.dedrs*drs2Ddxt) + d2facdt2*dtdxt*dtdxt*(-ggaw.f + ldaw.zk);

  r->d2fdxtxs[0] = 0.0;
  r->d2fdxtxs[1] = 0.0;
  r->d2fdxs2[0]  = 0.0;
  r->d2fdxs2[1]  = 0.0;
  r->d2fdxs2[2]  = 0.0;

}

#define func XC(gga_c_q2d_func)
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_q2d) = {
  XC_GGA_C_Q2D,
  XC_CORRELATION,
  "Chiodo et al",
  XC_FAMILY_GGA,
  "L Chiodo, LA Constantin, E Fabiano, and F Della Sala, Phys. Rev. Lett. 108, 126402 (2012)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-23, 0.0, 1e-32,
  gga_c_q2d_init,
  NULL, NULL,
  work_gga_c,
  NULL
};

