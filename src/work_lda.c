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
  This file is to be included in LDA functionals. As often these
  functionals are written as a function of rs and zeta, this
  routine performs the necessary conversions between this and a functional
  of rho.
************************************************************************/

static void 
work_lda(const void *p_, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2)
{
  const XC(lda_type) *p = p_;

  FLOAT dens, zeta, rs[3], drs;
  FLOAT zk0, dedrs, dedz, d2edrs2, d2edz2, d2edrsz;
  FLOAT *pdedrs, *pd2edrs2;

  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);

  /* Wigner radius */
  rs[1] = RS(dens);
  rs[0] = sqrt(rs[1]);
  rs[2] = rs[1]*rs[1];

  pdedrs   = (vrho!=NULL || v2rho2!=NULL) ? &dedrs   : NULL;
  pd2edrs2 =               (v2rho2!=NULL) ? &d2edrs2 : NULL;
  
  func(p, rs, zeta, &zk0, pdedrs, &dedz, pd2edrs2, &d2edrsz, &d2edz2);

  if(zk != NULL)
    *zk = zk0;

  if(vrho != NULL)
    vrho[0] = zk0 - (rs[1]/3.0)*dedrs;
    
  if(v2rho2 != NULL){
    drs = -(4.0*M_PI/9.0)*rs[2]*rs[2];
    v2rho2[0] = (2.0*dedrs - rs[1]*d2edrs2)*drs/3.0;
  }

  if(p->nspin == XC_POLARIZED){
    if(vrho != NULL){
      vrho[1] = vrho[0] - (zeta + 1.0)*dedz;
      vrho[0] = vrho[0] - (zeta - 1.0)*dedz;
    }

    if(v2rho2 != NULL){
      v2rho2[1] = v2rho2[0] - d2edrsz*(zeta - 1.0)*drs 
	+ (zeta + 1.0)/dens*(d2edrsz*rs[1]/3.0 + (zeta - 1.0)*d2edz2);

      v2rho2[2] = v2rho2[0] - d2edrsz*(zeta + 1.0)*drs
	+ (zeta + 1.0)/dens*(d2edrsz*rs[1]/3.0 + (zeta + 1.0)*d2edz2);
 
      v2rho2[0] = v2rho2[0] - d2edrsz*(zeta - 1.0)*drs 
	+ (zeta - 1.0)/dens*(d2edrsz*rs[1]/3.0 + (zeta - 1.0)*d2edz2);
      
    }
  }

}
