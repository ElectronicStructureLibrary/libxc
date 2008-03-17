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

#define XC_GGA_XC_XLYP 166 /* XLYP functional */

static void
gga_xc_xlyp_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  int i;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_x_init)(p->lda_aux, p->nspin, 3, XC_NON_RELATIVISTIC);

  p->gga_aux = (XC(gga_type) **) malloc(3*sizeof(XC(gga_type) *));
  for(i=0; i<3; i++)
    p->gga_aux[i] = (XC(gga_type) *) malloc(sizeof(XC(gga_type)));

  XC(gga_init)(p->gga_aux[0], XC_GGA_X_B88, p->nspin);
  XC(gga_init)(p->gga_aux[1], XC_GGA_X_PW91, p->nspin);
  XC(gga_init)(p->gga_aux[2], XC_GGA_C_LYP, p->nspin);
}

static void
gga_xc_xlyp_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  int i;

  for(i=0; i<3; i++){
    XC(gga_end)(p->gga_aux[i]);
    free(p->gga_aux[i]);
  }
  free(p->gga_aux);
}

static void 
gga_xc_xlyp(void *p_, FLOAT *rho, FLOAT *sigma,
            FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  static FLOAT cc[3] = {0.722, 0.347, 1.0};

  XC(gga_type) *p = p_;
  FLOAT dd, e1, vrho1[2], vsigma1[3];
  int ifunc, is, js;

  dd = 1.0 - cc[0] - cc[1];

  XC(lda_vxc)(p->lda_aux, rho, &e1, vrho1);
  *e = dd*e1;
  for(is=0; is<p->nspin; is++)
    vrho[is] = dd*vrho1[is];
  
  js = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
  for(is=0; is<js; is++)
    vsigma[is] = 0.0;

  for(ifunc=0; ifunc<3; ifunc++){
    XC(gga_vxc)(p->gga_aux[ifunc], rho, sigma, &e1, vrho1, vsigma1);

    *e += cc[ifunc]*e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += cc[ifunc]*vrho1[is];

    for(is=0; is<js; is++)
      vsigma[is] += cc[ifunc]*vsigma1[is];
  }
 
}


const XC(func_info_type) XC(func_info_gga_xc_xlyp) = {
  XC_GGA_XC_XLYP,
  XC_EXCHANGE_CORRELATION,
  "XLYP",
  XC_FAMILY_GGA,
  "X Xu and WA Goddard, III, PNAS 101, 2673 (2004)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_xlyp_init, 
  gga_xc_xlyp_end, 
  NULL,
  gga_xc_xlyp
};
