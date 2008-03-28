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

  FLOAT sfact, dens;
  int is;

  if(p->nspin == XC_POLARIZED){
    sfact     = 1.0;
  }else
    sfact     = 2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, ds, rho13;
    FLOAT x, f, dfdx, ldfdx, d2fdx2;
    FLOAT *pdfdx, *pvsigma, *pd2fdx2;
    int js = is==0 ? 0 : 2;

    if(rho[is] < MIN_DENS) continue;

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = POW(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);

    if(vrho!=NULL || v2rho2!=NULL){
      pdfdx   = &dfdx;
      pvsigma = &(vsigma[js]);
    }else
      pdfdx = pvsigma = NULL;

    d2fdx2 = 0.0; /* avoids a compiler warning */
    pd2fdx2 = (v2rho2!=NULL) ? &d2fdx2 : NULL;

#if   HEADER == 1
    func(p, x, &f, pdfdx, &ldfdx, pd2fdx2);
#elif HEADER == 2
    /* this second header is useful for functionals that depend
       explicitly both on s and on sigma */
    /* WARNING: This is not prepared for fxc */
    func(p, x, gdm*gdm, &f, pdfdx, &ldfdx, pvsigma);
#endif

    if(zk != NULL)
      *zk += -sfact*X_FACTOR_C*(ds*rho13)*f;
      
    if(vrho != NULL){
      vrho[is] += -4.0/3.0*X_FACTOR_C*rho13*(f - dfdx*x);
      if(gdm>MIN_GRAD)
	vsigma[js] = -sfact*X_FACTOR_C*(ds*rho13)*(vsigma[js]/(sfact*sfact) + dfdx*x/(2.0*sigma[js]));
      else
	vsigma[js] = -X_FACTOR_C/(sfact*(ds*rho13))*ldfdx;
    }

    if(v2rho2 != NULL){
      int n;

      v2rho2[js] = -4.0/3.0*X_FACTOR_C*rho13/(3.0*ds)*
	(f - dfdx*x + 4.0*d2fdx2*x*x)/sfact;

      n = (is == 0) ? 0 : 5;
      if(gdm>MIN_GRAD){
	v2rhosigma[n] =        X_FACTOR_C * 4.0/3.0*rho13 * d2fdx2*x*x/(2.0*sigma[js]);
	v2sigma2  [n] = -sfact*X_FACTOR_C*(ds*rho13)*
	  (d2fdx2*x - dfdx)*x/(4.0*sigma[js]*sigma[js]);
      }
	
    }
  }

  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}
