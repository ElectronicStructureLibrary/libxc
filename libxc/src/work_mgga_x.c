/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

/************************************************************************
  This file is to be included in meta GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3) and tau, this
  routine performs the necessary conversions between a functional of s and tau
  and of rho.
************************************************************************/

#include <stdio.h>

#ifndef XC_DIMENSIONS
#  define XC_DIMENSIONS 3
#endif

static void 
work_mgga_x(const void *p_, int np,
	    const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  FLOAT sfact, sfact2, dens, x_factor_c;
  int is, ip, order;
  int has_tail;

  /* WARNING: derivatives are _not_ OK for 2 dimensions */
  #if XC_DIMENSIONS == 2
  x_factor_c = X_FACTOR_2D_C;
  #else /* three dimensions */
  x_factor_c = X_FACTOR_C;
  #endif

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(order < 0) return;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  has_tail = 0;
  switch(p->info->number){
  case XC_MGGA_X_BR89:
  case XC_MGGA_X_BJ06:
  case XC_MGGA_X_TB09:
  case XC_MGGA_X_RPP09:
    has_tail = 1;
    break;
  }
  
  for(ip = 0; ip < np; ip++){
    dens = (p->nspin == XC_UNPOLARIZED) ? rho[0] : rho[0] + rho[1];
    if(dens < MIN_DENS) goto end_ip_loop;

    for(is=0; is<p->nspin; is++){
      FLOAT lrho, rho1D, rho2pD_D, lsigma, gdm;
      FLOAT x, t, u, f, lnr2, ltau, vrho0, dfdx, dfdt, dfdu;
      FLOAT d2fdx2, d2fdt2, d2fdu2, d2fdxt, d2fdxu, d2fdtu;
      int js = (is == 0) ? 0 : 2;
      int ks = (is == 0) ? 0 : 5;

      if((!has_tail && (rho[is] < MIN_DENS || tau[is] < MIN_TAU)) || (rho[is] == 0.0)) continue;

      lsigma= sigma[js]/sfact2;
      gdm   = sqrt(lsigma);
      lrho  = rho[is]/sfact;
      rho1D = POW(lrho, 1.0/XC_DIMENSIONS);
      rho2pD_D = lrho*rho1D*rho1D;
      x     = gdm/(lrho*rho1D);
    
      ltau  = tau[is]/sfact;
      t     = ltau/rho2pD_D;  /* tau/rho^((2+D)/D) */

      lnr2  = lapl_rho[is]/sfact;     /* this can be negative */
      u     = lnr2/rho2pD_D;  /* lapl_rho/rho^((2+D)/D) */

      vrho0 = dfdx = dfdt = dfdu = 0.0;
      d2fdx2 = d2fdt2 = d2fdu2 = d2fdxt = d2fdxu = d2fdtu = 0.0;

      func(p, x, t, u, order, &f, &vrho0,
	   &dfdx, &dfdt, &dfdu, 
	   &d2fdx2, &d2fdt2, &d2fdu2, &d2fdxt, &d2fdxu, &d2fdtu);

      if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
	*zk += -sfact*x_factor_c*(lrho*rho1D)*f;

      if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
	vrho[is]      = -x_factor_c*rho1D*(vrho0 + 4.0/3.0*(f - dfdx*x) - 5.0/3.0*(dfdt*t + dfdu*u));
	vtau[is]      = -x_factor_c*dfdt/rho1D;
	vlapl_rho[is] = -x_factor_c*dfdu/rho1D;
	if(gdm>MIN_GRAD)
	  vsigma[js] = -x_factor_c*(rho1D*lrho)*dfdx*x/(2.0*sfact*lsigma);
      }

      if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
	v2rho2[js]        = -x_factor_c/(9.0*sfact*rho1D*rho1D)*
	  (4.0*f - 4.0*x*dfdx + 4.0*4.0*x*x*d2fdx2 + 5.0*5.0*t*t*d2fdt2 + 5.0*5.0*u*u*d2fdu2 +
	   2.0*5.0*(4.0*x*t*d2fdxt + 4.0*x*u*d2fdxu + 5.0*t*u*d2fdtu));
	v2tau2[js]        = -x_factor_c*d2fdt2/(lrho*lrho);
	//v2lapl_rho2[js]   = -x_factor_c*d2fdu2/(rho1D*rho2pD_D);
	v2rhotau[js]      = -x_factor_c*rho1D/(3.0*sfact*rho2pD_D)*
	  (4.0*dfdt - 4.0*x*d2fdxt - 5.0*u*d2fdtu - 5.0*(dfdt + t*d2fdt2));
	//v2rholapl_rho[js] = -x_factor_c*rho1D/(3.0*sfact*rho2pD_D)*
	//  (4.0*dfdu - 4.0*x*d2fdxu - 5.0*u*d2fdtu - 5.0*(dfdu + u*d2fdu2));
	//v2taulapl_rho[js] = -x_factor_c*d2fdtu/(rho1D*rho2pD_D);
	if(gdm>MIN_GRAD){
	  v2sigma2[ks]   =  -x_factor_c*(rho1D*lrho)/(4.0*sfact2*sfact*lsigma*lsigma)*
	    (d2fdx2*x*x - dfdx*x);
	  v2rhosigma[ks] = -x_factor_c*rho1D*x/(3.0*2.0*sfact2*lsigma)*
	    (-4.0*x*d2fdx2 - 5.0*t*d2fdxt - 5.0*u*d2fdxu);
	  v2tausigma[ks] = -x_factor_c*x/(2.0*sfact2*lsigma*rho1D)*d2fdxt;
	  //v2sigmalapl_rho[ks] = -x_factor_c*x/(2.0*sfact2*lsigma*rho1D)*d2fdut;
	}
      }
    }
    
    if(zk != NULL)
      *zk /= dens; /* we want energy per particle */

  end_ip_loop:
    /* increment pointers */
    rho      += p->n_rho;
    sigma    += p->n_sigma;
    tau      += p->n_tau;
    lapl_rho += p->n_lapl_rho;
    
    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL){
      vrho      += p->n_vrho;
      vsigma    += p->n_vsigma;
      vtau      += p->n_vtau;
      vlapl_rho += p->n_vlapl_rho;
    }

    if(v2rho2 != NULL){
      v2rho2     += p->n_v2rho2;
      v2rhosigma += p->n_v2rhosigma;
      v2sigma2   += p->n_v2sigma2;
      /* warning: extra terms missing */
    }
  }
}
