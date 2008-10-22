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

/************************************************************************
  This file is to be included in meta GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3) and tau, this
  routine performs the necessary conversions between a functional of s and tau
  and of rho.
************************************************************************/

#include <stdio.h>

static void 
work_mgga_x(const void *p_, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
	    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  int is, order;

  order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, ds, rho13, sigmas;
    FLOAT x, t, f, ltau, dfdx, dfdt, d2fdx2, d2fdxt, d2fdt2;
    int js = (is == 0) ? 0 : 2;

    if(rho[is] < MIN_DENS) continue;

    dens  += rho[is];
    sigmas = max(MIN_GRAD*MIN_GRAD, sigma[js]/sfact2); 
    gdm    = sqrt(sigmas);
  
    ds    = rho[is]/sfact;
    rho13 = POW(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);
    
    ltau  = max(MIN_TAU, tau[is]/sfact);
    t     = ltau/(ds*rho13*rho13);  /* tau/rho^(5/3) */

    dfdx  = d2fdx2 = 0.0;

    dfdt = 0.0;
    func(p, x, t, order, &f, &dfdx, &dfdt, &d2fdx2, &d2fdxt, &d2fdt2);

    if(zk != NULL)
      *zk += -sfact*X_FACTOR_C*(ds*rho13)*f;

    if(vrho != NULL){
      vrho[is]   = -X_FACTOR_C*rho13*(4.0/3.0*(f - dfdx*x) - 5.0/3.0*dfdt*t);
      vtau[is]   = -X_FACTOR_C*dfdt/rho13;
      vsigma[js] = -X_FACTOR_C*(rho13*ds)*dfdx*x/(2.0*sfact*sigmas);
    }

    if(v2rho2 != NULL){
      /* Missing terms here */
      exit(1);
    }
  }

  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}
