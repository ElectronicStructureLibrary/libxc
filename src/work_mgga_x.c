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

static void 
work_mgga_x(const void *p_, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
	    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  int is;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, ds, rho13, z;
    FLOAT x, f, ltau, tauw, dfdx, dfdz, d2fdx2, d2fdxz, d2fdz2;
    FLOAT *pdfdx, *pd2fdx2;
    int js = (is == 0) ? 0 : 2;

    if(rho[is] < MIN_DENS) continue;

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = POW(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);
    
    ltau  = tau[is]/sfact;
    tauw  = gdm*gdm/(8.0*ds);
    z     = tauw/ltau;

    dfdx    = d2fdx2 = 0.0;
    pdfdx   = (vrho!=NULL || v2rho2!=NULL) ? &dfdx : NULL;
    pd2fdx2 = (v2rho2!=NULL) ? &d2fdx2 : NULL;

    dfdz = 0.0;
    func(p, x, z, &f, pdfdx, &dfdz, pd2fdx2, &d2fdxz, &d2fdz2);

    if(zk != NULL)
      *zk += -sfact*X_FACTOR_C*(ds*rho13)*f;
 
    if(vrho != NULL){
      vrho[is] += -X_FACTOR_C*rho13*(4.0/3.0*(f - dfdx*x) - dfdz*z);
      vtau[is] += -X_FACTOR_C*(ds*rho13)*dfdz*(-z/ltau);

      vsigma[js] = -sfact*X_FACTOR_C*rho13*
	(ds*dfdx*x/(2.0*sigma[js]) + dfdz/(8.0*ltau*sfact*sfact));
    }

    if(v2rho2 != NULL){
      /* Missing terms here */
      exit(1);
    }
  }

  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}
