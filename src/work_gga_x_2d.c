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

/************************************************************************
  This file is to be included in GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(3/2), this
  routine performs the necessary conversions between a functional of s
  and of rho.
************************************************************************/

/* TODO: merge this with 3D routine */

static void 
work_gga_x_2d(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const XC(gga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  int ip, is, order;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(order < 0) return;

  for(ip = 0; ip < np; ip++){
    dens = (p->nspin == XC_UNPOLARIZED) ? rho[0] : rho[0] + rho[1];
    if(dens < MIN_DENS) goto end_ip_loop;

    for(is=0; is<p->nspin; is++){
      FLOAT gdm, ds, sigmas, rho12;
      FLOAT x, f, dfdx, d2fdx2;
      int js = (is == 0) ? 0 : 2;
      int ks = (is == 0) ? 0 : 5;

      if(rho[is] < MIN_DENS) continue;

      sigmas = max(sigma[js], MIN_GRAD);
      gdm    = sqrt(sigmas)/sfact;

      ds    = rho[is]/sfact;
      rho12 = sqrt(ds);
      x     = gdm/(ds*rho12);

      f = dfdx = d2fdx2 = 0.0;

      func_2d(p, order, x, &f, &dfdx, &d2fdx2);

      if(zk != NULL)
	*zk += -sfact*X_FACTOR_2D_C*(ds*rho12)*f;
      
      if(vrho != NULL){
	vrho  [is] += -3.0/2.0*X_FACTOR_2D_C*rho12*(f - dfdx*x);
	vsigma[js]  = -sfact*X_FACTOR_2D_C*(ds*rho12)*dfdx*x/(2.0*sigmas);
      }
      
      if(v2rho2 != NULL){
	v2rho2[js] = -3.0/4.0*X_FACTOR_2D_C/rho12*
	  (f - dfdx*x + 3.0*d2fdx2*x*x)/sfact;
	
	if(gdm>MIN_GRAD){
	  v2rhosigma[ks] = -3.0/2.0*X_FACTOR_2D_C*rho12*(-d2fdx2*x*x/(2.0*sigmas));
	  v2sigma2  [ks] = -sfact*X_FACTOR_2D_C*(ds*rho12)*
	    (d2fdx2*x - dfdx)*x/(4.0*sigmas*sigmas);
	}
      }
    }
    
    if(zk != NULL)
      *zk /= dens; /* we want energy per particle */

  end_ip_loop:
    /* increment pointers */
    rho   += p->n_rho;
    sigma += p->n_sigma;
    
    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL){
      vrho   += p->n_vrho;
      vsigma += p->n_vsigma;
    }

    if(v2rho2 != NULL){
      v2rho2     += p->n_v2rho2;
      v2rhosigma += p->n_v2rhosigma;
      v2sigma2   += p->n_v2sigma2;
    }
  }
}
