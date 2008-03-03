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
#include <math.h>
#include "util.h"

void 
perdew_params(xc_gga_type *gga_p, FLOAT *rho, FLOAT *sigma, perdew_t *pt)
{
  pt->nspin = gga_p->nspin;
  rho2dzeta(pt->nspin, rho, &(pt->dens), &(pt->zeta));
  xc_lda_vxc(gga_p->lda_aux, rho, &(pt->ecunif), pt->vcunif);

  pt->rs = RS(pt->dens);
  pt->kf = POW(3.0*M_PI*M_PI*pt->dens, 1.0/3.0);
  pt->ks = sqrt(4.0*pt->kf/M_PI);

  pt->phi  = 0.5*(POW(1.0 + pt->zeta, 2.0/3.0) + POW(1.0 - pt->zeta, 2.0/3.0));

  /* get gdmt = |nabla n| */
  pt->gdmt = sigma[0];
  if(pt->nspin == XC_POLARIZED) pt->gdmt += 2.0*sigma[1] + sigma[2];
  if(pt->gdmt < MIN_GRAD*MIN_GRAD) pt->gdmt = MIN_GRAD*MIN_GRAD;
  pt->gdmt = sqrt(pt->gdmt);

  pt->t = pt->gdmt/(2.0 * pt->phi * pt->ks * pt->dens);

  pt->drs     = 0.0;
  pt->dkf     = 0.0;
  pt->dks     = 0.0;
  pt->dphi    = 0.0;
  pt->dt      = 0.0;
  pt->decunif = 0.0;
}

void 
perdew_potentials(perdew_t *pt, FLOAT *rho, FLOAT e_gga, 
		  FLOAT *vrho, FLOAT *vsigma)
{
  FLOAT drsdd, dkfdd, dksdd, dzdd[2], dpdz; 
  int is;
 
  drsdd   = -pt->rs/(3.0*pt->dens);
  dkfdd   =  pt->kf/(3.0*pt->dens);
  dksdd   = 0.5*pt->ks*dkfdd/pt->kf;
  dzdd[0] =  (1.0 - pt->zeta)/pt->dens;
  dzdd[1] = -(1.0 + pt->zeta)/pt->dens;
  dpdz    = 0.0;
  if(fabs(1.0 + pt->zeta) >= MIN_DENS)
    dpdz += (1.0/3.0)/POW(1.0 + pt->zeta, 1.0/3.0);
  if(fabs(1.0 - pt->zeta) >= MIN_DENS)
    dpdz -= (1.0/3.0)/POW(1.0 - pt->zeta, 1.0/3.0);
  
  /* add the t contributions to the other derivatives */
  pt->dphi += pt->dt * (-pt->t/pt->phi);
  pt->dks  += pt->dt * (-pt->t/pt->ks);

  /* calculate vrho */
  for(is=0; is<pt->nspin; is++){
    if(rho[is] > MIN_DENS){
      FLOAT decudd;
      
      vrho[is]  = e_gga;
      
      decudd = (pt->vcunif[is] - pt->ecunif)/pt->dens;
      vrho[is] += pt->dens * (  pt->decunif * decudd
			      + pt->drs     * drsdd
			      + pt->dkf     * dkfdd
			      + pt->dks     * dksdd
			      - pt->dt      * pt->t/pt->dens);
      vrho[is] += pt->dens * pt->dphi*dpdz*dzdd[is];
    }else{
      vrho[is] = 0.0;
    }
  }
    
  { /* calculate now vsigma */
    FLOAT dtdsig;
    
    dtdsig  = pt->t/(2.0*pt->gdmt*pt->gdmt);
    vsigma[0] = pt->dens*pt->dt*dtdsig;
    if(pt->nspin == XC_POLARIZED){
      vsigma[1] = 2.0*vsigma[0];
      vsigma[2] =     vsigma[0];
    }
  }
}
