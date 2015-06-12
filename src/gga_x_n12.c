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

#define XC_GGA_X_N12          82 /* N12 functional from Minnesota    */
#define XC_HYB_GGA_X_N12_SX   81 /* N12-SX functional from Minnesota */
#define XC_GGA_X_GAM          32 /* GAM functional from Minnesota */

static const FLOAT CC_N12[4][4] = {
  { 1.00000e+00,  5.07880e-01,  1.68233e-01,  1.28887e-01},
  { 8.60211e-02, -1.71008e+01,  6.50814e+01, -7.01726e+01},
  {-3.90755e-01,  5.13392e+01, -1.66220e+02,  1.42738e+02},
  { 4.03611e-01, -3.44631e+01,  7.61661e+01, -2.41834e+00}
};

static const FLOAT CC_N12_SX[4][4] = {
  /* Indices are wrong in the original paper; the first two indices
     need to be flipped */
  { 6.81116e-01,  1.88858e+00,  1.78590e+00,  8.79456e-01},
  {-8.12270e-02, -1.08723e+00, -4.18682e+00, -3.00000e+01},
  { 5.36236e-01, -5.45678e+00,  3.00000e+01,  5.51105e+01},
  {-7.09913e-01,  1.30001e+01, -7.24877e+01,  2.98363e+01}
};

static const FLOAT CC_GAM[4][4] = {
  { 1.32730,    0.886102, -5.73833,   8.60197},
  {-0.786018,  -4.78787,   3.90989,  -2.11611},
  { 0.802575,  14.4363,    8.842735, -6.21552},
  {-0.142331, -13.4598,    1.52355, -10.0530}
};

typedef struct{
  const FLOAT (*CC)[4];
} gga_x_n12_params;


static void
gga_x_n12_init(XC(func_type) *p)
{
  gga_x_n12_params *params;

  assert(p != NULL);

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_n12_params));
  params = (gga_x_n12_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_N12: 
    params->CC = CC_N12;
    break;
  case XC_HYB_GGA_X_N12_SX:
    params->CC = CC_N12_SX;
    p->cam_alpha = 0.00;
    p->cam_beta  = 0.25;
    p->cam_omega = 0.11;
    break;
  case XC_GGA_X_GAM:
    params->CC = CC_GAM;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_n12\n");
    exit(1);
  }
}


static void 
func(const XC(func_type) *pt, XC(gga_work_c_t) *r)
{
  gga_x_n12_params *params;

  int is;
  const FLOAT sign[2] = {1.0, -1.0}, omega_x=2.5, gamma_x=0.004;

  FLOAT opz, opz13, rss, x2;
  FLOAT vx, vx2, vx3, ux_d, ux, ux2, ux3;
  FLOAT pol1, pol2, pol3, pol4;
  FLOAT ex, FN12;

  FLOAT drssdrs, drssdz, dvxdrss, duxdxs;
  FLOAT dpol1, dpol2, dpol3, dpol4;
  FLOAT dexdz, dexdrss, dFN12dux, dFN12dvx;

  assert(pt != NULL && pt->params != NULL);
  params = (gga_x_n12_params *) (pt->params);

  r->f = 0.0;
  if(r->order >= 1)
    r->dfdrs = r->dfdz = r->dfdxt = r->dfdxs[0] = r->dfdxs[1] = 0.0;

  /* now the spin-resolved part */
  for(is = 0; is < 2; is++){
    opz   = 1.0 + sign[is]*r->zeta;
    if(opz < pt->info->min_zeta) continue;

    opz13 = CBRT(opz);
    rss   = r->rs*M_CBRT2/opz13;
    x2    = r->xs[is]*r->xs[is];

    vx    = 1.0/(1.0 + (1.0/(RS_FACTOR*omega_x))*rss);

    ux_d  = 1.0/(1.0 + gamma_x*x2);
    ux    = gamma_x*x2*ux_d;

    vx2 = vx*vx; vx3 = vx2*vx;
    ux2 = ux*ux; ux3 = ux2*ux;

    pol1 = params->CC[0][0] + params->CC[0][1]*ux + params->CC[0][2]*ux2 + params->CC[0][3]*ux3;
    pol2 = params->CC[1][0] + params->CC[1][1]*ux + params->CC[1][2]*ux2 + params->CC[1][3]*ux3;
    pol3 = params->CC[2][0] + params->CC[2][1]*ux + params->CC[2][2]*ux2 + params->CC[2][3]*ux3;
    pol4 = params->CC[3][0] + params->CC[3][1]*ux + params->CC[3][2]*ux2 + params->CC[3][3]*ux3;

    FN12 = pol1 + vx*pol2 + vx2*pol3 + vx3*pol4;

    ex    = -X_FACTOR_C*RS_FACTOR*opz/(2.0*rss);
    r->f += ex*FN12;

    if(r->order < 1) continue;

    drssdrs = M_CBRT2/opz13;
    drssdz  = -sign[is]*rss/(3.0*opz);

    dvxdrss = -vx*vx/(RS_FACTOR*omega_x);
    duxdxs  = 2.0*gamma_x*r->xs[is]*ux_d*ux_d;

    dpol1 = params->CC[0][1] + 2.0*params->CC[0][2]*ux + 3.0*params->CC[0][3]*ux2;
    dpol2 = params->CC[1][1] + 2.0*params->CC[1][2]*ux + 3.0*params->CC[1][3]*ux2;
    dpol3 = params->CC[2][1] + 2.0*params->CC[2][2]*ux + 3.0*params->CC[2][3]*ux2;
    dpol4 = params->CC[3][1] + 2.0*params->CC[3][2]*ux + 3.0*params->CC[3][3]*ux2;

    dFN12dux = dpol1 + vx*dpol2 + vx2*dpol3 + vx3*dpol4;
    dFN12dvx = pol2 + 2.0*vx*pol3 + 3.0*vx2*pol4;

    dexdrss = -ex/rss;
    dexdz   = sign[is]*ex/opz;

    r->dfdrs    += (dexdrss*FN12 + ex*dFN12dvx*dvxdrss)*drssdrs;
    r->dfdz     += dexdz*FN12 + (dexdrss*FN12 + ex*dFN12dvx*dvxdrss)*drssdz;
    r->dfdxs[is] = ex*dFN12dux*duxdxs;
  }
}


#include "work_gga_c.c"


const XC(func_info_type) XC(func_info_gga_x_n12) = {
  XC_GGA_X_N12,
  XC_EXCHANGE,
  "Minnesota N12 exchange functional to be used with gga_c_n12",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2012_2310, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_x_n12_init,
  NULL, NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_x_n12_sx) = {
  XC_HYB_GGA_X_N12_SX,
  XC_EXCHANGE,
  "Worker for hyb_gga_x_n12_sx",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Peverati2012_16187, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_x_n12_init,
  NULL, NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_gam) = {
  XC_GGA_X_GAM,
  XC_EXCHANGE,
  "GAM functional from Minnesota",
  XC_FAMILY_GGA,
  {&xc_ref_Yu2015_12146, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_x_n12_init,
  NULL, NULL,
  work_gga_c,
  NULL
};
