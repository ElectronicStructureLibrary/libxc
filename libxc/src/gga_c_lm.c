/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
  This functional is provided for historical reasons.
  It was one of the first GGAs that ever appeared.
************************************************************************/

#define XC_GGA_C_LM          137 /* Langreth and Mehl correlation          */

static void 
gga_c_lm_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_init)(p->lda_aux, XC_LDA_C_vBH, p->nspin);
}


static void 
gga_c_lm_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  free(p->lda_aux);
}


static void 
gga_c_lm(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	 FLOAT *e, FLOAT *vrho, FLOAT *vsigma,
	 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  XC(perdew_t) pt;

  int order;
  FLOAT me;
  FLOAT H, dHdx1, d2Hdx12, dx1dt, dx1dphi;
  FLOAT grad_to_t;
  FLOAT x1;

  const FLOAT a1 = 4.28e-3/2.0; /* The factor of 2 converts from Rydberg to Hartree */
  const FLOAT a2 = -0.262;
  const FLOAT a3 = 7.0/9.0;

  order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  XC(perdew_params)(p, rho, sigma, order, &pt);

  grad_to_t = sqrt(4.0/M_PI)*POW(3*M_PI*M_PI, 1.0/6.0);
  x1 = 2.0*grad_to_t*pt.phi*pt.t;

  H = a1*x1*x1*(exp(a2*x1) - a3);

  me = pt.ecunif + H;
  if(e != NULL) *e = me;

  if(order >= 1){
    dHdx1      = a1*x1*((2.0 + x1*a2)*exp(a2*x1) - 2.0*a3);
    dx1dphi    = 2.0*grad_to_t*  pt.t;
    dx1dt      = 2.0*grad_to_t*pt.phi;

    pt.dphi    = dHdx1 * dx1dphi;
    pt.dt      = dHdx1 * dx1dt;
    pt.decunif = 1.0;
  }

  if(order >= 2){
    d2Hdx12 = a1*((2.0 + a2*x1*(4.0 + a2*x1))*exp(a2*x1) - 2.0*a3);

    pt.d2phi2 = d2Hdx12*dx1dphi*dx1dphi;
    pt.d2phit = 2.0*grad_to_t*dHdx1 + d2Hdx12*dx1dphi*dx1dt;
    pt.d2t2   = d2Hdx12*dx1dt*dx1dt;
  }

  XC(perdew_potentials)(&pt, rho, me, order, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}


const XC(func_info_type) XC(func_info_gga_c_lm) = {
  XC_GGA_C_LM,
  XC_CORRELATION,
  "Langreth & Mehl",
  XC_FAMILY_GGA,
  "DC Langreth and MJ Mehl, Phys. Rev. Lett. 47, 446 (1981)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  gga_c_lm_init,
  gga_c_lm_end,
  NULL,            /* this is not an LDA                   */
  gga_c_lm,
};

