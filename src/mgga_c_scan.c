/*
 Copyright (C) 2016 Susi Lehtola

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

#define XC_MGGA_C_SCAN          267 /* SCAN correlation */

/* Constants */
const static FLOAT b1c = 0.0285764, b2c = 0.0889, b3c = 0.125541;

static void mgga_c_scan_init(XC(func_type) *p)
{
  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_C_SCAN_E0,    p->nspin);
}

/* Calculates E_c^{LDA0} = -b1c / ( 1 + b2c*r^{1/2} + b3c*r ) */
static void
func_eclda0(FLOAT rs, int order, FLOAT *e, FLOAT *dedrs)
{
  FLOAT sqrtrs, den, dden;

  sqrtrs = SQRT(rs);
  den  = 1.0 + b2c*sqrtrs + b3c*rs;
  
  *e = -b1c / den;

  if(order < 1) return;

  dden   = b2c/(2.0*sqrtrs) + b3c;
  *dedrs = DFRACTION(-b1c, 0.0, den, dden);
}

/* Calculates g(chi=infty,s) = 1 / ( 1 + 4 chi(infty) s^2 )^{1/4} */
static void
func_g_infty(FLOAT s, int order, FLOAT *g, FLOAT *dgds)
{
  const static FLOAT chi_infty = 0.12802585262625815;
  FLOAT aux;

  aux = 1 + 4.0*chi_infty*s*s;
  *g = 1.0 / SQRT(SQRT(aux));

  if(order < 1) return;

  *dgds = -2.0*chi_infty*s*(*g)/aux;
}

/* Calculates d and phi: [ (1+z)^a + (1-z)^a ] / 2 */
/* This function should go global and replace eventually FZETA */
static void
func_dnphi(FLOAT z, FLOAT a, int order, FLOAT *phi, FLOAT *dphidz)
{
  FLOAT opza, omza;

  opza = POW(1.0+z,a);
  omza = POW(1.0-z,a);

  *phi = 0.5*(opza + omza);

  if(order < 1) return;

  *dphidz = ABS(z)==1.0 ? (FLT_MAX) : 0.5*a*(opza/(1.0 + z) - omza/(1.0 - z));
}


/* Calculates Gc = [ 1 - 2.3621(dx(z) - 1)] (1 - z^12) */
static void
func_Gc(FLOAT z, int order, FLOAT *Gc, FLOAT *dGcdz)
{
  const static FLOAT G_c = 2.363; /* in the paper it is 2.3631 */
  FLOAT dx, ddxdz, zp11;
  
  zp11 = POW(z, 11);
  func_dnphi(z, 4.0/3.0, order, &dx, &ddxdz);

  *Gc = (1.0 - G_c*(dx - 1.0))*(1.0 - z*zp11);

  if(order < 1) return;

  *dGcdz = -G_c*ddxdz*(1.0 - z*zp11) - 12.0*zp11*(1.0 - G_c*(dx - 1.0));
}

/* Calculates 
   H0 = b1c ln [ 1 + w0(rs) (1 - g(s)) ]
*/
static void
func_e0(FLOAT rs, FLOAT zeta, FLOAT s, int order, 
        FLOAT *e0, FLOAT *de0drs, FLOAT *de0dz, FLOAT *de0ds)
{
  FLOAT expn, eclda0, declda0drs, w0, dw0drs, dw0dz;
  FLOAT gs, dgds, Gc, dGcdz, H0aux, H0, dH0drs, dH0ds;
  
  func_eclda0(rs, order, &eclda0, &declda0drs);
  
  expn = EXP(-eclda0/b1c);
  w0   = expn - 1;

  func_g_infty(s, order, &gs, &dgds);
  
  H0aux = 1.0 + w0*(1.0 - gs);
  H0   = b1c*LOG(H0aux);

  func_Gc(zeta, order, &Gc, &dGcdz);

  *e0 = (eclda0 + H0)*Gc;

  if(order < 1) return;

  dw0drs = -declda0drs*expn/b1c;

  dH0drs =  b1c*(1.0 - gs)*dw0drs/H0aux;
  dH0ds  = -b1c*w0*dgds/H0aux;

  *de0drs = (declda0drs + dH0drs)*Gc;
  *de0dz  = (eclda0 + H0)*dGcdz;
  *de0ds  = dH0ds*Gc;
}


static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  const static FLOAT c1c=0.64, c2c=1.5, dc=0.7;

  XC(gga_work_c_t) pbe;
  FLOAT ss, e0, de0drs, de0dz, de0ds;
  FLOAT opz, omz, opz13, omz13, opz23, omz23;
  FLOAT a_num, a_den, da_numdz, da_dendz, alpha, dadz, dadxt, dadts[2], fc, dfcda;

  pbe.order = r->order;
  pbe.rs    = r->rs;
  pbe.z     = r->zeta;
  pbe.xt    = r->xt;
  pbe.xs[0] = r->xs[0]; 
  pbe.xs[1] = r->xs[1];

  XC(gga_c_scan_e0_func) (pt->func_aux[0], &pbe);

  ss = X2S*M_CBRT2*r->xt;
  func_e0(r->rs, r->zeta, ss, r->order, &e0, &de0drs, &de0dz, &de0ds);

  opz = 1.0 + r->zeta;
  omz = 1.0 - r->zeta;

  opz13 = CBRT(opz);
  omz13 = CBRT(omz);

  opz23 = opz13*opz13;
  omz23 = omz13*omz13;

  a_num = r->ts[0]*opz*opz23 + r->ts[1]*omz*omz23 - r->xt*r->xt*M_CBRT2*M_CBRT2/4.0;
  a_den = K_FACTOR_C*(opz*opz23 + omz*omz23);
  alpha = a_num / a_den;

  XC(mgga_x_scan_falpha)(r->order, alpha, c1c, c2c, dc, &fc, &dfcda);

  r->f = pbe.f + fc*(e0 - pbe.f);

  if(r->order < 1) return;

  da_numdz = (5.0/3.0)*(r->ts[0]*opz23 - r->ts[1]*omz23);
  da_dendz = (K_FACTOR_C*5.0/3.0)*(opz23 - omz23);

  dadz     = DFRACTION(a_num, da_numdz, a_den, da_dendz);
  dadxt    = -r->xt*M_CBRT2*M_CBRT2/(2.0*a_den);
  dadts[0] = opz*opz23/a_den;
  dadts[1] = omz*omz23/a_den;

  r->dfdrs = pbe.dfdrs + fc*(de0drs - pbe.dfdrs);
  r->dfdz  = pbe.dfdz  + fc*(de0dz  - pbe.dfdz) + (e0 - pbe.f)*dfcda*dadz;
  r->dfdxt = pbe.dfdxt + fc*(de0ds*X2S*M_CBRT2 - pbe.dfdxt) + (e0 - pbe.f)*dfcda*dadxt;
  r->dfdxs[0] = r->dfdxs[1] = 0.0;
  r->dfdts[0] = (e0 - pbe.f)*dfcda*dadts[0];
  r->dfdts[1] = (e0 - pbe.f)*dfcda*dadts[1];
}


#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_scan) = {
  XC_MGGA_C_SCAN,
  XC_CORRELATION,
  "SCAN correlation of Sun, Ruzsinszky, and Perdew",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_DEVELOPMENT,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  mgga_c_scan_init, 
  NULL, NULL, NULL,
  work_mgga_c,
};
