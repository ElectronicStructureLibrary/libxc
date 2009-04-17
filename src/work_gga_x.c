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
  functionals are written as a function of s = |grad n|/n^(4/3), this
  routine performs the necessary conversions between a functional of s
  and of rho.
************************************************************************/

#ifndef HEADER
#  define HEADER 1
#endif

static void 
work_gga_x(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const XC(gga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  int is, order;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(order < 0) return;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, ds, rho13;
    FLOAT x, f, dfdx, ldfdx, d2fdx2, lvsigma, lv2sigma2, lvsigmax;
    int js = (is == 0) ? 0 : 2;
    int ks = (is == 0) ? 0 : 5;

    if(rho[is] < MIN_DENS) continue;

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = POW(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);

    dfdx = ldfdx = d2fdx2 = 0.0;
#if   HEADER == 1
    func(p, order, x, &f, &dfdx, &ldfdx, &d2fdx2);
    lvsigma = lv2sigma2 = lvsigmax = 0.0;
#elif HEADER == 2
    /* this second header is useful for functionals that depend
       explicitly both on s and on sigma */
    lvsigma = lv2sigma2 = lvsigmax = 0.0;
    func(p, order, x, gdm*gdm, &f, &dfdx, &ldfdx, &lvsigma, &d2fdx2, &lv2sigma2, &lvsigmax);
#endif

    lvsigma   /= sfact2;
    lvsigmax  /= sfact2;
    lv2sigma2 /= sfact2*sfact2;

    if(zk != NULL && (p->info->provides & XC_PROVIDES_EXC))
      *zk += -sfact*X_FACTOR_C*(ds*rho13)*f;
      
    if(vrho != NULL && (p->info->provides & XC_PROVIDES_VXC)){
      vrho[is] += -4.0/3.0*X_FACTOR_C*rho13*(f - dfdx*x);
      if(gdm>MIN_GRAD)
	vsigma[js] = -sfact*X_FACTOR_C*(ds*rho13)*(lvsigma + dfdx*x/(2.0*sigma[js]));
      else
	vsigma[js] = -X_FACTOR_C/(sfact*(ds*rho13))*ldfdx;
    }

    if(v2rho2 != NULL && (p->info->provides & XC_PROVIDES_FXC)){
      v2rho2[js] = -4.0/3.0*X_FACTOR_C*rho13/(3.0*ds)*
	(f - dfdx*x + 4.0*d2fdx2*x*x)/sfact;

      if(gdm>MIN_GRAD){
	v2rhosigma[ks] = -4.0/3.0*X_FACTOR_C*rho13*(lvsigma - lvsigmax*x - d2fdx2*x*x/(2.0*sigma[js]));
	v2sigma2  [ks] = -sfact*X_FACTOR_C*(ds*rho13)*
	  (lv2sigma2 + lvsigmax*x/sigma[js] + (d2fdx2*x - dfdx)*x/(4.0*sigma[js]*sigma[js]));
      }
	
    }
  }

  if(zk != NULL && (p->info->provides & XC_PROVIDES_EXC))
    *zk /= dens; /* we want energy per particle */
}
