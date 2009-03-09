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
work_lda(const void *p_, const FLOAT *rho, 
	 FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
{
  const XC(lda_type) *p = p_;

  XC(lda_rs_zeta) r;
  int is;
  FLOAT dens, dcnst, drs, rs4, d2rs;

  r.order = -1;
  if(zk     != NULL) r.order = 0;
  if(vrho   != NULL) r.order = 1;
  if(v2rho2 != NULL) r.order = 2;
  if(v3rho3 != NULL) r.order = 3;
  if(r.order < 0) return;

  XC(rho2dzeta)(p->nspin, rho, &dens, &r.zeta);

  /* Wigner radius */
  r.rs[1] = RS(dens);
  r.rs[0] = sqrt(r.rs[1]);
  r.rs[2] = r.rs[1]*r.rs[1];

  func(p, &r);

  if(zk != NULL && (p->info->provides & XC_PROVIDES_EXC))
    *zk = r.zk;

  if(r.order < 1) return;

  if(vrho != NULL && (p->info->provides & XC_PROVIDES_VXC)){
    vrho[0] = r.zk - (r.rs[1]/3.0)*r.dedrs;

    if(p->nspin == XC_POLARIZED){
      vrho[1] = vrho[0] - (r.zeta + 1.0)*r.dedz;
      vrho[0] = vrho[0] - (r.zeta - 1.0)*r.dedz;
    }
  }
  
  if(r.order < 2) return;

  dcnst = 4.0*M_PI/3.0;
  rs4   = r.rs[2]*r.rs[2];
  drs   = -dcnst*rs4/3.0;

  if(v2rho2 != NULL && (p->info->provides & XC_PROVIDES_FXC)){
    v2rho2[0] = (2.0*r.dedrs - r.rs[1]*r.d2edrs2)*drs/3.0;
    
    if(p->nspin == XC_POLARIZED){
      FLOAT sign[3][2] = {{-1.0, -1.0}, {-1.0, +1.0}, {+1.0, +1.0}};

      for(is=2; is>=0; is--)
	v2rho2[is] = v2rho2[0] - r.d2edrsz*(r.zeta + sign[is][0])*drs
	  + (r.zeta + sign[is][1])/dens*(r.d2edrsz*r.rs[1]/3.0 + (r.zeta + sign[is][0])*r.d2edz2);
    }
  }

  if(r.order < 3) return;

  if(v3rho3 != NULL && (p->info->provides & XC_PROVIDES_KXC)){
    d2rs = (4.0/9.0) * dcnst*dcnst * rs4*r.rs[2]*r.rs[1];

    v3rho3[0]  = (r.d2edrs2 - r.rs[1]*r.d3edrs3)*drs*drs + (2.0*r.dedrs - r.rs[1]*r.d2edrs2)*d2rs;
    v3rho3[0] /= 3.0;

    if(p->nspin == XC_POLARIZED){
      FLOAT sign[4][3] = {{-1.0, -1.0, -1.0}, {-1.0, -1.0, +1.0}, {-1.0, +1.0, +1.0}, {+1.0, +1.0, +1.0}};

      for(is=3; is>=0; is--){
	FLOAT ff;

	v3rho3[is]  = v3rho3[0] - (r.zeta + sign[is][0])*(r.d3edrs2z*drs*drs + r.d2edrsz*d2rs);

	v3rho3[is] += (r.zeta + sign[is][1])/dens*
	  (r.d3edrs2z*r.rs[1]/3.0 + r.d2edrsz/3.0 + (r.zeta + sign[is][0])*r.d3edrsz2)*drs;

	v3rho3[is] += -(r.zeta + sign[is][1])/(dens*dens)*
	  (r.d2edrsz*r.rs[1]/3.0 + (r.zeta + sign[is][0])*r.d2edz2);

	ff  = (2.0*r.d2edrsz - r.rs[1]*r.d3edrs2z)*drs/3.0;
	ff += -(r.d3edrsz2*(r.zeta + sign[is][0]) + r.d2edrsz)*drs;
	ff += (r.d2edrsz*r.rs[1]/3.0 + (r.zeta + sign[is][0])*r.d2edz2)/dens;
	ff += (r.zeta + sign[is][1])/dens*(r.d3edrsz2*r.rs[1]/3.0 + r.d2edz2 + (r.zeta + sign[is][0])*r.d3edz3);

	v3rho3[is] += -ff*(r.zeta + sign[is][2])/dens;
      }
    }
  }
}
