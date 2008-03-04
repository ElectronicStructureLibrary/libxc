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

#include "util.h"

#define XC_MGGA_C_TPSS          202 /* Perdew, Tao, Staroverov & Scuseria correlation */

/* WARNING - this is all broken !! */
#define   _(is, x)   [3*is + x]
#define  __(i, j)    [2*i + j] 


/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.
   J. Chem. Phys. 120, 6898 (2004)
   http://dx.doi.org/10.1063/1.1665298

  Correlation part
************************************************************************/

static XC(func_info_type) func_info_mgga_c_tpss = {
  XC_MGGA_C_TPSS,
  XC_CORRELATION,
  "Perdew, Tao, Staroverov & Scuseria",
  XC_FAMILY_MGGA,
  "J.P.Perdew, Tao, Staroverov, and Scuseria, Phys. Rev. Lett. 91, 146401 (2003)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC
};


void XC(mgga_c_tpss_init)(XC(mgga_type) *p)
{
  p->info = &func_info_mgga_c_tpss;

  p->gga_aux1 = (XC(gga_type) *) malloc(sizeof(XC(gga_type)));
  XC(gga_init)(p->gga_aux1, XC_GGA_C_PBE, p->nspin);

  if(p->nspin == XC_UNPOLARIZED){
    p->gga_aux2 = (XC(gga_type) *) malloc(sizeof(XC(gga_type)));
    XC(gga_init)(p->gga_aux2, XC_GGA_C_PBE, XC_POLARIZED);
  }
}


void XC(mgga_c_tpss_end)(XC(mgga_type) *p)
{
  XC(gga_end)(p->gga_aux1);
  free(p->gga_aux1);

  if(p->nspin == XC_UNPOLARIZED) {
    XC(gga_end)(p->gga_aux2);
    free(p->gga_aux2);
  }
}


/* some parameters */
static FLOAT d = 2.8;


/* Equation (14) */
static void
c_tpss_14(FLOAT csi, FLOAT zeta, FLOAT *C, FLOAT *dCdcsi, FLOAT *dCdzeta)
{
  FLOAT fz, C0, dC0dz, dfzdz;
  FLOAT z2 = zeta*zeta;
    
  /* Equation (13) */
  C0    = 0.53 + z2*(0.87 + z2*(0.50 + z2*2.26));
  dC0dz = zeta*(2.0*0.87 + z2*(4.0*0.5 + z2*6.0*2.26));
  
  fz    = 0.5*(POW(1.0 + zeta, -4.0/3.0) + POW(1.0 - zeta, -4.0/3.0));
  dfzdz = 0.5*(POW(1.0 + zeta, -7.0/3.0) - POW(1.0 - zeta, -7.0/3.0))*(-4.0/3.0);
  
  { /* Equation (14) */
    FLOAT csi2 = csi*csi;
    FLOAT a = 1.0 + csi2*fz, a4 = POW(a, 4);
    
    *C      =  C0 / a4;
    *dCdcsi = -8.0*csi*fz/(a*a4);
    *dCdzeta = (dC0dz*a - C0*4.0*csi2*dfzdz)/(a*a4);
  }
}


/* Equation 12 */
static void c_tpss_12(XC(mgga_type) *p, FLOAT *rho, FLOAT *grho, 
		 FLOAT dens, FLOAT zeta, FLOAT z,
		 FLOAT *e_PKZB, FLOAT *de_PKZBdd, FLOAT *de_PKZBdgd, FLOAT *de_PKZBdz)
{
  FLOAT e_PBE, *de_PBEdd, *de_PBEdgd;
  FLOAT e_til[2], de_tildd[2], de_tildgd[2*3];

  FLOAT C, dCdcsi, dCdzeta;
  FLOAT *dzetadd, *dcsidd, *dcsidgd;
  int i, is;

  de_PBEdd  = (FLOAT *)malloc(p->nspin*sizeof(FLOAT));
  de_PBEdgd = (FLOAT *)malloc(3*p->nspin*sizeof(FLOAT));
  dzetadd   = (FLOAT *)malloc(p->nspin*sizeof(FLOAT));
  dcsidd    = (FLOAT *)malloc(p->nspin*sizeof(FLOAT));
  dcsidgd   = (FLOAT *)malloc(3*p->nspin*sizeof(FLOAT));

  { /* get the PBE stuff */
    XC(gga_type) *aux2 = (p->nspin == XC_UNPOLARIZED) ? p->gga_aux2 : p->gga_aux1;
    XC(gga)(p->gga_aux1, rho, grho, &e_PBE, de_PBEdd, de_PBEdgd);
    
    for(is=0; is<p->nspin; is++){
      FLOAT r1[2], gr1[2*3], e1, de1dd[2], de1dgd[2*3];
      FLOAT sfac = (p->nspin == XC_UNPOLARIZED) ? 0.5 : 1.0;
      
      /* build fully polarized density and gradient */
      r1[0] = rho[is] * sfac;
      r1[1] = 0.0;
      
      for(i=0; i<3; i++){
	gr1 _(0, i) = grho _(is, i) * sfac;
	gr1 _(1, i) = 0.0;
      }

      /* call (polarized) PBE again */
      de1dd[0] = 0.0; de1dd[1] = 0.0;
      XC(gga)(aux2, r1, gr1, &e1, de1dd, de1dgd);

      e_til   [is] = e1;
      de_tildd[is] = de1dd[0];
      for(i=0; i<3; i++) de_tildgd _(is, i) = de1dgd _(0, i);
    }
  } /* end PBE stuff */


  if(p->nspin == XC_UNPOLARIZED){
    C          = 0.53;
    dzetadd[0] = 0.0;
    dcsidd [0] = 0.0;
    for(i=0; i<3; i++) dcsidgd _(0, i) = 0.0;

  }else{ /* get C(csi, zeta) */
    FLOAT gzeta[3], gzeta2, csi, a;
    
    for(i=0, gzeta2=0.0; i<3; i++){
      gzeta[i] = 2.0*(grho _(0, i)*rho[1] - grho _(1, i)*rho[0])/(dens*dens);
      gzeta2  += gzeta[i]*gzeta[i];
    }
    gzeta2 = max(gzeta2, MIN_GRAD*MIN_GRAD);
    
    a = 2.0*POW(3.0*M_PI*M_PI*dens, 1.0/3.0);
    csi = sqrt(gzeta2)/a;
  
    c_tpss_14(csi, zeta, &C, &dCdcsi, &dCdzeta);
    
    dzetadd[0] =  (1.0 - zeta)/dens;
    dzetadd[1] = -(1.0 + zeta)/dens;

    dcsidd [0] = -7.0*csi/(3.0*dens); 
    dcsidd [1] = dcsidd [0];
    for(i=0; i<3; i++){
      FLOAT aa = gzeta[i]/(sqrt(gzeta2)*a);
      dcsidd[0] += -grho _(1, i)*aa;
      dcsidd[1] +=  grho _(0, i)*aa;
      
      dcsidgd _(0, i) =  rho[1]*aa;
      dcsidgd _(1, i) = -rho[0]*aa;
    }    
  } /* get C(csi, zeta) */

  { /* end */
    FLOAT z2 = z*z, aux, *dauxdd, *dauxdgd;

    dauxdd  = (FLOAT *)malloc(p->nspin*sizeof(FLOAT));
    dauxdgd = (FLOAT *)malloc(3*p->nspin*sizeof(FLOAT));

    /* aux = sum_sigma n_sigma e_til */
    aux = 0.0;
    for(is=0; is<p->nspin; is++){
      aux += rho[is] * max(e_til[is], e_PBE);

      dauxdd[is] = 0.0;
    }

    for(is=0; is<p->nspin; is++){
      int is2 = (is == 0) ? 1 : 0;
      
      dauxdd[is] -= aux/dens;
      
      if(e_til[is] > e_PBE){
	dauxdd[is]  += e_til[is] + rho[is]*de_tildd[is]/dens;
	for(i=0; i<3; i++){
	  dauxdgd _(is, i) += rho[is]*de_tildgd _(is, i)/dens;
	}
      }else{
	dauxdd[is]  += e_PBE + rho[is] * de_PBEdd[is];
	dauxdd[is2] +=         rho[is] * de_PBEdd[is2];
	for(i=0; i<3; i++){
	  dauxdgd _(is, i)  += rho[is]*de_PBEdgd _(is,  i)/dens;
	  dauxdgd _(is2, i) += rho[is]*de_PBEdgd _(is2, i)/dens;
	}
      }
    }
 
    *e_PKZB    = e_PBE*(1 + C*z2) - (1.0 + C)*z2*aux;
    *de_PKZBdz = dens * 2.0*z * C * (e_PBE - aux);
    for(is=0; is<p->nspin; is++){
      FLOAT dCdd;
      
      dCdd = dCdzeta*dzetadd[is] + dCdcsi*dcsidd[is];
      
      de_PKZBdd[is] = de_PBEdd[is]*(1.0 + C*z2) + dens*e_PBE*dCdd*z2;
      de_PKZBdd[is]-= z2*(dCdd*aux + (1.0 + C)*dauxdd[is]);
			  
      for(i=0; i<3; i++){
	FLOAT dCdgd =  dCdcsi*dcsidgd _(is, i);
	
	de_PKZBdgd _(is, i) = de_PBEdgd _(is, i)*(1.0 + C*z2) + dens*e_PBE*dCdgd*z2;
	de_PKZBdgd _(is, i)-= z2*(dCdgd*aux + (1.0 + C)*dauxdgd _(is, i));
	
      }
    }
  } /* end */
}


void 
XC(mgga_c_tpss)(XC(mgga_type) *p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
	    FLOAT *energy, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau)
{
  FLOAT dens, zeta;
  FLOAT gd[3], gdms, taut, tauw, z;
  FLOAT e_PKZB, *de_PKZBdd, *de_PKZBdgd, de_PKZBdz;
  int i, is;

  de_PKZBdd  = (FLOAT *)malloc(p->nspin*sizeof(FLOAT));
  de_PKZBdgd = (FLOAT *)malloc(3*p->nspin*sizeof(FLOAT));

  /* change variables */
  XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);

  /* sum tau and grho over spin */
  for(i=0; i<3; i++) gd[i] = grho _(0, i);
  taut = tau[0];

  for(is=1; is<p->nspin; is++){
    for(i=0; i<3; i++) gd[i] += grho _(is, i);
    taut += tau[is];
  }

  /* get the modulos square of the gradient */
  gdms = gd[0]*gd[0] + gd[1]*gd[1] + gd[2]*gd[2];
  gdms = max(MIN_GRAD*MIN_GRAD, gdms);
  
  tauw = gdms/(8.0*dens);

  /* sometimes numerical errors makes taut negative :( */
  if(taut < MIN_TAU)
    taut = tauw;

  z = tauw/taut;

  /* Equation (12) */
  c_tpss_12(p, rho, grho, z, dens, zeta,
	    &e_PKZB, de_PKZBdd, de_PKZBdgd, &de_PKZBdz);
  
  /* Equation (11) */
  {
    FLOAT z2 = z*z, z3 = z2*z;
    FLOAT dedz;
    FLOAT dzdd, dzdgd[3], dzdtau;

    dzdd   = -z/dens;
    dzdtau = -z/taut;
    for(i=0; i<3; i++) dzdgd[i] = gd[i]/(4.0*dens*taut);

    *energy = e_PKZB*(1.0 + d*e_PKZB*z3);
    dedz    = de_PKZBdz*(1.0 + 2.0*d*e_PKZB*z3) +
      dens * e_PKZB*e_PKZB * d * 3.0*z2;  

    for(is=0; is<p->nspin; is++){
      dedd[is]   = de_PKZBdd[is] * (1.0 + 2.0*d*e_PKZB*z3);
      dedd[is]  -= e_PKZB*e_PKZB * d * z3;
      dedd[is]  += dedz*dzdd;
      
      for(i=0; i<3; i++){
	dedgd _(is, i) = de_PKZBdgd _(is, i) * (1.0 + 2.0*d*e_PKZB*z3);
	dedgd _(is, i)+= dedz*dzdgd[i];
      }
      
      dedtau[is] = dedz*dzdtau;
    }
  }

}
