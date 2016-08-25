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
#define XC_HYB_MGGA_X_SCAN0     264 /* SCAN hybrid exchange */

static void
func_gx(int order, FLOAT ss, FLOAT *gx, FLOAT *dgxds)
{
  static const FLOAT a1  = 4.9479;
  FLOAT thr, smh, smhps, expn;

  thr = a1*a1/(LOG(FLOAT_EPSILON)*LOG(FLOAT_EPSILON));

  /* Special handling for small values of s */
  if(ss < thr) {
    smh = smhps = expn = 0.0;
  } else {
    smh   = 1.0/SQRT(ss);
    smhps = smh/ss;
    expn  = EXP(-a1*smh);
  }

  *gx = 1.0 - expn;

  if(order < 1) return;

  *dgxds = -0.5 * a1 * expn * smhps;
}

void
XC(mgga_x_scan_falpha)(int order, FLOAT a, FLOAT c1, FLOAT c2, FLOAT dd, FLOAT *f, FLOAT *dfda)
{
  /* exponentials are truncated */
  const FLOAT logeps =  LOG(FLOAT_EPSILON);
  FLOAT thr1, thr2;
  FLOAT c1exp, c2exp, ooma;

  thr1  = -logeps/(c1 - logeps);
  thr2  = 1.0 - c2/logeps;

  ooma = 1.0/(1.0 - a);

  c1exp = (a >= thr1) ? 0.0 : EXP(-c1*a*ooma);
  c2exp = (a <= thr2) ? 0.0 : EXP(c2*ooma);

  *f = c1exp - dd*c2exp;

  if(order < 1) return;

  *dfda = -(c1*c1exp + dd*c2*c2exp)*ooma*ooma;
}

static void
func_x(int order, FLOAT ss, FLOAT a, FLOAT *x, FLOAT *dxds, FLOAT *dxda)
{
  static const FLOAT mu  = 10.0/81.0;
  static const FLOAT b1  = 0.1566320774354852; /* (511.0/13500.0)/(2*b2) */
  static const FLOAT b2  = 0.1208304597359457; /* SQRT(5913.0/405000.0) */
  static const FLOAT b3  = 0.5;
  static const FLOAT b4  = 0.1218315102059958; /* mu^2/k1 - 1606.0/18225.0 - b1^2 */

  FLOAT ss2, beta, expb3, expb4, sterm;

  ss2  = ss*ss;
  beta = 1.0 - a;

  /* Helpers */
  expb4 = EXP(-b4*ss2/mu);
  expb3 = EXP(-b3*beta*beta);

  /* Second term in the bracket */
  sterm = b1*ss2 + b2*beta*expb3;

  *x = mu*ss2 + b4*ss2*ss2*expb4 + sterm*sterm;

  if(order < 1) return;

  *dxds = 2.0*ss*(mu + b4*ss2*(2.0 - b4*ss2/mu)*expb4 + 2*b1*(b1*ss2 + b2*beta*expb3));
  *dxda = -2.0*(b1*ss2 + b2*beta*expb3) * b2 * expb3 * (1.0 - 2.0*b3*beta*beta);
}

static void
func_h1x(int order, FLOAT x, FLOAT *h, FLOAT *dhdx)
{
  static const FLOAT k1  = 0.065;
  FLOAT k1ok1px;

  k1ok1px = k1/(k1 + x);
  *h = 1.0 + k1*(1.0 - k1ok1px);

  if(order < 1) return;

  *dhdx = k1ok1px*k1ok1px;
}

static void
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  static const FLOAT c1x = 0.667, c2x = 0.8, dx  = 1.24;
  static const FLOAT h0x = 1.174;

  /* s variable and alpha */
  FLOAT s, a;
  /* Derivatives of alpha in terms of libxc variables */
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

  s = X2S*r->x;
  a = (r->t - r->x*r->x/8.0)/K_FACTOR_C;

  /* Calculate functions */
  XC(mgga_x_scan_falpha)(r->order, a, c1x, c2x, dx, &fx, &dfxda);
  func_gx(r->order, s, &gx, &dgxds);
  func_x(r->order, s, a, &y, &dyds, &dyda);
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
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  NULL, NULL, NULL, NULL,
  work_mgga_x,
};

static void
hyb_mgga_x_scan0_init(XC(func_type) *p)
{
  static int   funcs_id  [1] = {XC_MGGA_X_SCAN};
  static FLOAT funcs_coef[1] = {1.0 - 0.25};

  XC(mix_init)(p, 1, funcs_id, funcs_coef);
  p->cam_alpha = 0.25;
}


const XC(func_info_type) XC(func_info_hyb_mgga_x_scan0) = {
  XC_HYB_MGGA_X_SCAN0,
  XC_EXCHANGE,
  "SCAN hybrid exchange (SCAN0)",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Hui2016_044114, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  hyb_mgga_x_scan0_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};
