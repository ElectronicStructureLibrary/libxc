/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_X_TAU_HCTH        205 /* tau-HCTH from Boese and Handy */
#define XC_MGGA_X_BMK             279 /* Boese-Martin for kinetics     */
#define XC_HYB_MGGA_X_TAU_HCTH    282 /* Hybrid version of tau-HCTH    */

const FLOAT tHCTH_cx_local [4] = {1.10734, -1.0534, 6.3491, -2.5531};
const FLOAT tHCTH_cx_nlocal[4] = {0.00110, -0.3041, 6.9543, -0.7235};

const FLOAT BMK_cx_local [4] = { 0.474302, 2.77701, -11.4230, 11.7167};
const FLOAT BMK_cx_nlocal[4] = {-0.192212, 4.73936, -26.6188, 22.4891};

const FLOAT hyb_tHCTH_cx_local [4] = { 0.86735,  0.3008, 1.2208,   0.1574};
const FLOAT hyb_tHCTH_cx_nlocal[4] = {-0.00230, -0.2849, 5.4146, -10.909};

typedef struct{
  const FLOAT *cx_local;
  const FLOAT *cx_nlocal;
} mgga_x_tau_hcth_params;


static void 
mgga_x_tau_hcth_init(XC(func_type) *p)
{
  mgga_x_tau_hcth_params *params;

  assert(p != NULL);
  assert(p->params == NULL);

  p->params = malloc(sizeof(mgga_x_tau_hcth_params));
  params = (mgga_x_tau_hcth_params *)(p->params);

  switch(p->info->number){
  case XC_MGGA_X_TAU_HCTH:
    params->cx_local  = tHCTH_cx_local;
    params->cx_nlocal = tHCTH_cx_nlocal;
    break;
  case XC_MGGA_X_BMK:
    params->cx_local  = BMK_cx_local;
    params->cx_nlocal = BMK_cx_nlocal;
    break;
  case XC_HYB_MGGA_X_TAU_HCTH:
    p->cam_alpha = 0.15;
    params->cx_local  = hyb_tHCTH_cx_local;
    params->cx_nlocal = hyb_tHCTH_cx_nlocal;
    break;
  default:
    fprintf(stderr, "Internal error in mgga_tau_hcth\n");
    exit(1);
    break;
  }
}


/* Eq. (22) */
void
XC(mgga_b00_fw)(int order, FLOAT t, FLOAT *fw, FLOAT *dfwdt)
{
  FLOAT w, w2;

  w = (K_FACTOR_C - t)/(K_FACTOR_C + t);
  w2 = w*w;

  *fw = w*(1.0 - 2.0*w2 + w2*w2);

  if(order < 1) return;

  *dfwdt = 1.0 - 6.0*w2 + 5.0*w2*w2;
  *dfwdt *= -2.0*K_FACTOR_C/((K_FACTOR_C + t)*(K_FACTOR_C + t));
}


static void
eq_29(int order, FLOAT x, FLOAT *ux, FLOAT *duxdx)
{
  static FLOAT gamX = 0.004;
  FLOAT x2, denom;

  x2    = x*x;
  denom = 1.0 + gamX*x2;

  *ux = gamX*x2/denom;

  if(order < 1) return;

  *duxdx = 2.0*gamX*x/(denom*denom);
}

static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  const mgga_x_tau_hcth_params *params;
  const FLOAT *cx_local, *cx_nlocal;

  FLOAT ux, ux2, gxl, gxnl, fx;
  FLOAT duxdx, dgxldu, dgxnldu, dfxdt;

  params = (mgga_x_tau_hcth_params *)(pt->params);
  cx_local  = params->cx_local;
  cx_nlocal = params->cx_nlocal;

  eq_29          (r->order, r->x, &ux, &duxdx);
  XC(mgga_b00_fw)(r->order, r->t, &fx, &dfxdt);

  ux2  = ux*ux;
  gxl  = cx_local [0] + ux*(cx_local [1] + cx_local [2]*ux + cx_local [3]*ux2);
  gxnl = cx_nlocal[0] + ux*(cx_nlocal[1] + cx_nlocal[2]*ux + cx_nlocal[3]*ux2);

  r->f = gxl + gxnl*fx;

  if(r->order < 1) return;

  dgxldu  = cx_local [1] + 2.0*cx_local [2]*ux + 3.0*cx_local [3]*ux2;
  dgxnldu = cx_nlocal[1] + 2.0*cx_nlocal[2]*ux + 3.0*cx_nlocal[3]*ux2;

  r->dfdx = (dgxldu + dgxnldu*fx)*duxdx;
  r->dfdt = gxnl*dfxdt;
}

#include "work_mgga_x.c"

const XC(func_info_type) XC(func_info_mgga_x_tau_hcth) = {
  XC_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "tau-HCTH from Boese and Handy",
  XC_FAMILY_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  mgga_x_tau_hcth_init, 
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_mgga_x_bmk) = {
  XC_MGGA_X_BMK,
  XC_EXCHANGE,
  "Boese-Martin for kinetics",
  XC_FAMILY_MGGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  mgga_x_tau_hcth_init, 
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_hyb_mgga_x_tau_hcth) = {
  XC_HYB_MGGA_X_TAU_HCTH,
  XC_EXCHANGE,
  "Hybrid version of tau-HCTH",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Boese2002_9559, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  mgga_x_tau_hcth_init, 
  NULL, NULL, NULL,
  work_mgga_x,
};
