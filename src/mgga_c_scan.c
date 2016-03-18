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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_MGGA_C_SCAN          267 /* SCAN correlation */

/* Constants */
const FLOAT b1c=0.02858, b2c=0.0889, b3c=0.1255;
const FLOAT c1c=0.64, c2c=1.5, dc=0.7;
const FLOAT c_gamma=0.03109069086965490; /* (1-log(2))/Pi^2 */
const FLOAT beta_a=0.066725, beta_b=0.1, beta_c=0.1778;
const FLOAT G_c=2.3621;
const FLOAT chi=0.128026;

/* Calculates E_c^{LDA0} = -b1c / ( 1 + b2c*r^{1/2} + b3c*r ) */
static void
Eclda0(FLOAT r, int order,
       FLOAT *e, FLOAT *dedr)
{
  FLOAT sqrtr=SQRT(r);
  FLOAT denom=1.0 + b2c*sqrtr + b3c*r;
  
  *e = -b1c / denom;

  if(order < 1) return;

  *dedr = *e/denom * (b2c/(2.0*sqrtr) + b3c);
}

/* Calculates g(x) = 1 / ( 1 + 4 c x^2 )^{1/4} */
static void
func_gx(FLOAT x, FLOAT c, int order,
       FLOAT *g, FLOAT *dgdx)
{
  *g = 1.0 / SQRT(SQRT(1 + 4.0*c*x*x));

  if(order < 1) return;

  *dgdx = -2.0*c*x/(1 + 4.0*c*x*x)*(*g);
}

/* Calculates d and phi: [ (1+z)^a + (1-z)^a ] / 2 */
static void
func_dnphi(FLOAT z, FLOAT a, int order,
	FLOAT *phi, FLOAT *dphidz)
{
  FLOAT opza=POW(1.0+z,a);
  FLOAT omza=POW(1.0-z,a);

  *phi = 0.5*(opza + omza);

  if(order < 1) return;

  *dphidz = 0.5*a*(opza/(1.0+z) - omza/(1.0-z));
}

/* Calculate d_x(z) */
static void
func_d(FLOAT z, int order,
	FLOAT *d, FLOAT *dddz)
{
  func_dnphi(z,4.0/3.0,order,d,dddz);
}

/* Calculate phi(z) */
static void
func_phi(FLOAT z, int order,
	FLOAT *phi, FLOAT *dphidz)
{
  func_dnphi(z,2.0/3.0,order,phi,dphidz);
}

/* Calculates Gc = [ 1 - 2.3621( d_x(z)-1 )] (1 - z^12) */
static void
func_Gc(FLOAT z, int order,
	FLOAT *G, FLOAT *dGdz)
{
  /* d_x and dg_x/dz */
  FLOAT dx, ddxdz;
  /* z^11 */
  FLOAT zp11 = POW(z,11);

  /* Calculate d_x */
  func_phi(z,order,&dx,&ddxdz);

  *G = (1.0 - G_c*(dx - 1.0))*(1.0 - z*zp11);

  if(order < 1) return;

  *dGdz = -G_c*ddxdz*(1.0 - z*zp11) - 12.0*zp11*(1.0 - G_c*(dx - 1.0));
}

/* Calculates w(r,z) = exp( - E_c^{LDA0}(r) / f(z) ) */
static void
func_w(FLOAT r, int order,
       FLOAT fz, FLOAT dfdz,
       FLOAT *w, FLOAT *dwdr, FLOAT *dwdz)
{
  FLOAT expn;
  FLOAT E0, dE0dr;
  
  /* Get E_c^{LDA0} */
  Eclda0(r,order,&E0,&dE0dr);

  expn = EXP(-E0/fz);
  *w = expn - 1;

  if(order < 1) return;

  *dwdr = -dE0dr/fz*expn;
  *dwdz = E0*dfdz/(fz*fz)*expn;
}

/* Calculates beta(r) = 0.066725 (1 + 0.1 r_s) / (1 + 0.1778 r_s) */
static void
func_beta(FLOAT r, int order,
       FLOAT *b, FLOAT *dbdr)
{
  FLOAT denom=1.0 + beta_c*r;
  
  *b = beta_a * (1.0 + beta_b*r) / denom;

  if(order < 1) return;

  *dbdr = beta_a * (beta_b - beta_c) / (denom*denom);
}

/* Calculates A(r,f(z)) = beta(r)/f(z) */
static void
func_A(FLOAT r, int order,
       FLOAT fz, FLOAT dfdr, FLOAT dfdz,
       FLOAT *A, FLOAT *dAdr, FLOAT *dAdz)
{
  /* Calculate beta */
  FLOAT b, dbdr;
  func_beta(r,order,&b,&dbdr);
  
  *A = b/fz;

  if(order < 1) return;

  *dAdr = dbdr/fz - b*dfdr/(fz*fz);
  *dAdz = -b*dfdz/(fz*fz);
}

/* Calculates t = (3 pi^2 / 16)^{1/3} s / (phi r^1/2) */
static void
func_t(FLOAT r, FLOAT s, int order,
       FLOAT phi, FLOAT dphidz,
       FLOAT *t, FLOAT *dtdr, FLOAT *dtdz, FLOAT *dtds)
{
  FLOAT c=CBRT(3.0*M_PI*M_PI/16.0);
  FLOAT rh=SQRT(r);
  
  *t = c*s / (phi*rh);
  
  if(order < 1) return;

  *dtds = c / (phi*rh);
  *dtdr = -(*t)/(2*r);
  *dtdz = -(*t)*dphidz/phi;
}

/* Calculates 
   H0 = b1c ln [ 1 + w0(r) (1 - g(s)) ]
*/
static void
func_H0(FLOAT r, FLOAT s,
	int order, int n,
	FLOAT *H, FLOAT *dHdr, FLOAT *dHds)
{
  /* w_0 and its derivatives */
  FLOAT w, dwdr, dwdz;
  /* g and its derivatives */
  FLOAT gs, dgds;
  /* derivative denominator */
  FLOAT ddenom;
  
  /* Calculate w */
  func_w(r,order,b1c,0.0,&w,&dwdr,&dwdz);

  /* Calculate g */
  func_gx(s,chi,order,&gs,&dgds);
  
  /* Value of H is */
  *H = b1c * LOG( 1.0 + w * (1 - gs));

  if(order < 1) return;

  /* Derivatives */
  ddenom = b1c/(1.0 + w*(1-gs));
  *dHds = ddenom*(-w*dgds);
  *dHdr = ddenom*dwdr*(1-gs);
}

/* Calculates 
   H1 = (gamma phi^3) ln [ 1 + w1(r,z) (1 - g(At^2)) ]
*/
static void
func_H1(FLOAT r, FLOAT z, FLOAT s,
	int order, int n,
	FLOAT *H, FLOAT *dHdr, FLOAT *dHdz, FLOAT *dHds)
{
  /* Phi and its derivative */
  FLOAT phi, dphidz;
  FLOAT phi2, phi3;
  /* w_1 and its derivatives */
  FLOAT argw, dargwdz;
  FLOAT w, dwdr, dwdz;
  /* A and its derivatives */
  FLOAT argA, dargAdr, dargAdz;
  FLOAT A, dAdr, dAdz;
  /* t and its derivatives */
  FLOAT t, dtdr, dtdz, dtds;
  /* g and its derivatives */
  FLOAT x, dxdr, dxdz, dxds;
  FLOAT gx, dgdx;
  FLOAT dgdr, dgdz, dgds;
  /* logarithmic term */
  FLOAT ln;
  /* derivative denominator */
  FLOAT ddenom;
  
  /* Calculate phi */
  func_phi(z,order,&phi,&dphidz);
  phi2=phi*phi;
  phi3=phi2*phi;
  
  /* Calculate w. Argument is */
  argw=c_gamma*phi3;
  /* and its derivative */
  dargwdz=3*c_gamma*phi2*dphidz;
  /* and w is */
  func_w(r,order,argw,dargwdz,&w,&dwdr,&dwdz);

  /* Calculate A */
  argA=c_gamma*w;
  dargAdr=c_gamma*dwdr;
  dargAdz=c_gamma*dwdz;
  func_A(r,order,argA,dargAdr,dargAdz,&A,&dAdr,&dAdz);

  /* Calculate t */
  func_t(r,s,order,phi,dphidz,&t,&dtdr,&dtdz,&dtds);

  /* Calculate g */
  x = A*t*t;
  func_gx(x,1.0,order,&gx,&dgdx);
  
  /* Value of H is */
  ln = LOG( 1.0 + w * (1 - gx));
  *H = argw * ln;

  if(order < 1) return;

  /* Derivatives of g */
  dgdr = dgdx*(2.0*A*t*dtdr);
  dgdz = dgdx*(2.0*A*t*dtdz);
  dgds = dgdx*(2.0*A*t*dtds);
  
  /* Derivatives */
  ddenom = argw/(1.0 + w*(1-gx));
  *dHds = ddenom*(-w*dgds);
  *dHdz = dargwdz*ln + ddenom*(dwdz*(1-gx) - w*dgdz);
  *dHdr = ddenom*(dwdr*(1-gx) - w*dgdr);
}

static void
func_fx(int order, FLOAT a, FLOAT *f, FLOAT *dfda)
{
  FLOAT c1exp=0.0, c2exp=0.0;
  FLOAT ooma=1.0/(1.0-a);

  c1exp=XC(mgga_x_scan_exp1)(c1c,a);
  c2exp=XC(mgga_x_scan_exp2)(c2c,a);
  *f = c1exp - dc*c2exp;

  if(order < 1) return;

  *dfda = -(c1c*c1exp + dc*c2c*c2exp)*ooma*ooma;
}

static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  XC(lda_work_t) pw;
  pw.order = r->order;
  pw.rs[0] = SQRT(r->rs);
  pw.rs[1] = r->rs;
  pw.rs[2] = r->rs*r->rs;
  pw.zeta  = r->zeta;

  XC(lda_c_pw_func)(pt->func_aux[0], &pw);

  fprintf(stderr, "SCAN correlation not implemented yet.\n");
  exit(1);
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
  NULL, NULL, NULL, NULL,
  work_mgga_c,
};
