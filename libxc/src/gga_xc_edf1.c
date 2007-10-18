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

#define XC_GGA_XC_EDF1 165 /* Empirical functionals from Adamson, Gill, and Pople */

void gga_xc_edf1_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;
  int i;

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_x_init(p->lda_aux, p->nspin, 3, XC_NON_RELATIVISTIC);

  p->gga_aux = (xc_gga_type **) malloc(3*sizeof(xc_gga_type *));
  for(i=0; i<3; i++)
    p->gga_aux[i] = (xc_gga_type *) malloc(sizeof(xc_gga_type));

  xc_gga_init(p->gga_aux[0], XC_GGA_X_B88, p->nspin);
  gga_x_b88_set_params(p->gga_aux[0], 0.0035);

  xc_gga_init(p->gga_aux[1], XC_GGA_X_B88, p->nspin);
  gga_x_b88_set_params(p->gga_aux[1], 0.0042);

  xc_gga_init(p->gga_aux[2], XC_GGA_C_LYP, p->nspin);
  gga_c_lyp_set_params(p->gga_aux[2], 0.055, 0.158, 0.25, 0.3505);
}

void gga_xc_edf1_end(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;
  int i;

  for(i=0; i<3; i++){
    xc_gga_end(p->gga_aux[i]);
    free(p->gga_aux[i]);
  }
  free(p->gga_aux);
}

static void 
gga_xc_edf1(void *p_, double *rho, double *sigma,
            double *e, double *vrho, double *vsigma)
{
  static double cx    = 1.030952;
  static double cc[3] = {10.4017, -8.44793, 1.0};

  xc_gga_type *p = p_;
  double dd, e1, vrho1[2], vsigma1[3];
  int ifunc, is, js;

  xc_lda_vxc(p->lda_aux, rho, &e1, vrho1);
  dd = cx - cc[0] - cc[1];
  *e = dd*e1;
  for(is=0; is<p->nspin; is++)
    vrho[is] = dd*vrho1[is];
  
  js = (p->nspin == XC_UNPOLARIZED) ? 1 : 3;
  for(is=0; is<js; is++)
    vsigma[is] = 0.0;

  for(ifunc=0; ifunc<3; ifunc++){
    xc_gga(p->gga_aux[ifunc], rho, sigma, &e1, vrho1, vsigma1);

    *e += cc[ifunc]*e1;
    for(is=0; is<p->nspin; is++)
      vrho[is] += cc[ifunc]*vrho1[is];

    //if(ifunc == 2) continue;
    for(is=0; is<js; is++)
      vsigma[is] += cc[ifunc]*vsigma1[is];
  }
  
}

const xc_func_info_type func_info_gga_xc_edf1 = {
  XC_GGA_XC_EDF1,
  XC_EXCHANGE_CORRELATION,
  "EDF1",
  XC_FAMILY_GGA,
  "RD Adamson, PMW Gill, and JA Pople, Chem. Phys. Lett. 284 6 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_edf1_init, 
  gga_xc_edf1_end, 
  NULL,
  gga_xc_edf1
};
