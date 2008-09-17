/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

#include <stdlib.h>
#include <assert.h>
#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_LTA          201 /* Local tau approximation of Ernzerhof & Scuseria */

static void
mgga_x_lta(const void *p_, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
	   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau,
	   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  /*C_x* POW(10.0/(3.0*POW(3.0*M_PI*M_PI, 2.0/3.0)), 4.0/5.0) */
  const FLOAT a1 = -X_FACTOR_C*0.430075922439080216009;

  FLOAT sfact, ltau, dens;
  int is;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    int js = (is == 0) ? 0 : 2;

    dens += rho[is];
    ltau  = tau[is]/sfact;

    if(zk != NULL)
      *zk += a1*POW(ltau, 4.0/5.0);

    if(vrho != NULL)
      vtau[is] = a1*4.0/5.0*POW(ltau, -1.0/5.0);

    if(v2rho2 != NULL)
      v2tau2[js] = -a1*4.0/25.0*POW(ltau, -6.0/5.0)/sfact;
  }
  
  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}

const XC(func_info_type) XC(func_info_mgga_x_lta) = {
  XC_MGGA_X_LTA,
  XC_EXCHANGE,
  "Local tau approximation",
  XC_FAMILY_MGGA,
  "M Ernzerhof and G Scuseria, J. Chem. Phys. 111, 911 (1999)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL, NULL,
  NULL, NULL,        /* this is not an LDA                   */
  mgga_x_lta,
};
