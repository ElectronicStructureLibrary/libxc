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

#define XC_MGGA_X_SCAN          263 /* SCAN exchange of Sun, Ruzsinszky, and Perdew  */


static void
func_gx(int order, FLOAT s, FLOAT *g, FLOAT *dgds)
{
  const FLOAT a1=4.9479;
  FLOAT smh, smhps, expn;

  /* Special handling for small values of s */
  const FLOAT thr=a1*a1/(LOG(FLOAT_EPSILON)*LOG(FLOAT_EPSILON));

  if(s < thr) {
    smh=0.0;
    smhps=0.0;
    expn=0.0;
  } else {
    smh=1.0/SQRT(s);
    smhps=smh/s;
    expn=EXP(-a1*smh);
  }
    
  *g = 1.0 - expn;

  if(order < 1) return;

  *dgds = - 0.5 * a1 * expn * smhps;
}

static void
func_fx(int order, FLOAT a, FLOAT *f, FLOAT *dfda)
{
  const FLOAT c1x=0.667;
  const FLOAT c2x=0.8;
  const FLOAT dx=1.24;

  FLOAT c1exp=0.0, c2exp=0.0;
  FLOAT dc1exp=0.0, dc2exp=0.0;
  
  /* Step functions implemented here */
  if(a < 1.0)
    c1exp=EXP(-c1x*a/(1.0-a));
  else if(a > 1.0)
    c2exp=EXP(-c2x/(1.0-a));
  
  *f = c1exp - dx*c2exp;

  if(order < 1) return;

  /* Step functions implemented here */
  if(a < 1.0)
    dc1exp=-c1x*c1exp/((1.0-a)*(1.0-a));
  else if(a > 1.0)
    dc2exp=-c2x*c1exp/((1.0-a)*(1.0-a));
  
  *dfda = dc1exp - dx*dc2exp;
}

  
static void
func_x(int order, FLOAT s, FLOAT a,
       FLOAT *x, FLOAT *dxds, FLOAT *dxda)
{
  const FLOAT mu=10.0/81.0;
  const FLOAT b2=SQRT(5913.0/405000.0);
  const FLOAT b1=(511.0/13500.0)/(2*b2);
  const FLOAT b3=0.5;
  const FLOAT k1=0.065;
  const FLOAT b4=mu*mu/k1 - 1606.0/18225.0 - b1*b1;

  /* Write in terms of variables */
  FLOAT p=s*s;
  FLOAT beta=1.0-a;

  /* Helpers */
  FLOAT expb4=EXP(-b4*p/mu);
  FLOAT expb3=EXP(-b3*beta*beta);
  /* Second term in the bracket */
  FLOAT sterm=b1*p + b2*beta*expb3;

  *x = mu*p + b4*p*p*expb4 + sterm*sterm;

  if(order < 1) return;

  *dxds = 2.0*s*(mu + b4*p*(2.0 - b4*p/mu)*expb4 + 2*b1*(b1*p + b2*beta*expb3));
  *dxda = -2.0*( b1*p + b2*beta*expb3 ) * b2 * expb3 * (1.0 - 2.0*b3*beta*beta);
}

static void
func_h1x(int order, FLOAT x,
	 FLOAT *h, FLOAT *dhdx)
{
  const FLOAT k1=0.065;

  *h = 1.0 + k1*(1.0 - k1/(k1+x));
  
  if(order < 1) return;

  *dhdx = k1*k1/((k1+x)*(k1+x));
}

static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  /* s variable and alpha */
  FLOAT ss, a;
    /* Derivatives of alpha */
  FLOAT dadx, dadt;

  /* x is an internal variable in libxc which differs from the x used in the functional.
     For clarity, denote x of the functional with y. y and its derivatives */
  FLOAT y, dyds, dyda;
  
  /* h1x(y) and its derivative */
  FLOAT h1x, dh1xdy;

  /* f_x(alpha) */
  FLOAT fx, dfxda;

  /* g_x(s) */
  FLOAT gx, dgxds;

  /* Derivatives of full functional */
  FLOAT dFds, dFda;
  
  /* h0x */
  const FLOAT h0x=1.174;

  ss = X2S*r->x;
  a = (r->t - r->x*r->x/8.0)/K_FACTOR_C;

  /* Calculate functions */
  func_fx(r->order, a, &fx, &dfxda); 
  func_gx(r->order, ss, &gx, &dgxds); 
  func_x(r->order, ss, a, &y, &dyds, &dyda);
  func_h1x(r->order, y, &h1x, &dh1xdy); 
  
  /* Functional value is */
  r->f = (h1x + fx*(h0x-h1x))*gx;

  if(r->order < 1) return;

  dFds = (1.0 - fx)*gx*dh1xdy*dyds + (h1x + fx*(h0x-h1x))*dgxds;
  dFda = ((1.0 - fx)*dh1xdy*dyda + dfxda*(h0x-h1x))*gx;
  dadx = -2.0*r->x/(8.0*K_FACTOR_C);
  dadt = 1.0/K_FACTOR_C;

  r->dfdx = dFda*dadx + dFds*X2S;
  r->dfdt = dFda*dadt;
  r->dfdu = 0.0;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_scan) = {
  XC_MGGA_X_SCAN,
  XC_EXCHANGE,
  "SCAN exchange of Sun, Ruzsinszky, and Perdew",
  XC_FAMILY_MGGA,
  {&xc_ref_Sun2015_036402, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_DEVELOPMENT,
  1e-32, 1e-32, 1e-32, 1e-32,
  NULL, NULL, NULL, NULL,
  work_mgga_x,
};
