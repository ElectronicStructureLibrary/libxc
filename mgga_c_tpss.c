#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Perdew, Tao, Staroverov & Scuseria 
   meta-Generalized Gradient Approximation.

  Correlation part
************************************************************************/

static func_type func_mgga_c_tpss = {
  XC_MGGA_C_TPSS,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof",
  "MGGA",
  {"J.P.Perdew, Tao, Staroverov, and Scuseria, Phys. Rev. Lett. 91, 146401 (2003)", NULL}
};


void mgga_c_tpss_init(mgga_type *p)
{
  p->func = &func_mgga_c_tpss;

  p->gga_aux1 = (gga_type *) malloc(sizeof(gga_type));
  gga_init(p->gga_aux1, XC_GGA_C_PBE, p->nspin);

  if(p->nspin == XC_UNPOLARIZED){
    p->gga_aux2 = (gga_type *) malloc(sizeof(gga_type));
    gga_init(p->gga_aux1, XC_GGA_C_PBE, XC_POLARIZED);
  }
}


void mgga_c_tpss_end(mgga_type *p)
{
  gga_end(p->gga_aux1);
  free(p->gga_aux1);

  if(p->nspin == XC_UNPOLARIZED) {
    gga_end(p->gga_aux2);
    free(p->gga_aux2);
  }
}


/* some parameters */
static double d=2.8;


/* Equation (14) */
static void
c_tpss_C(double csi, double zeta, double *C, double *dCdcsi, double *dCdzeta)
{
  double fz, C0, dC0dz, dfzdz;
  double z2 = zeta*zeta;
    
  /* Equation (13) */
  C0    = 0.53 + z2*(0.87 + z2*(0.50 + z2*2.26));
  dC0dz = zeta*(2.0*0.87 + z2*(4.0*0.5 + z2*6.0*2.26));
  
  fz    = 0.5*(pow(1.0 + zeta, -4.0/3.0) + pow(1.0 - zeta, -4.0/3.0));
  dfzdz = 0.5*(pow(1.0 + zeta, -7.0/3.0) - pow(1.0 - zeta, -7.0/3.0))*(-4.0/3.0);
  
  { /* Equation (14) */
    double csi2 = csi*csi;
    double a = 1.0 + csi2*fz, a4 = pow(a, 4);
    
    *C      =  C0 / a4;
    *dCdcsi = -8.0*csi*fz/(a*a4);
    *dCdzeta = (dC0dz*a - C0*4.0*csi2*dfzdz)/(a*a4);
  }
}


/* Equation 12 */
static void c_tpss_12(mgga_type *p, double *rho, double *grho, 
		 double dens, double zeta, double z,
		 double *e_PKZB, double *de_PKZBdd, double *de_PKZBdgd, double *de_PKZBdz)
{
  double e_PBE,           de_PBEdd[p->nspin], de_PBEdgd[p->nspin*3];
  double e_til[p->nspin], de_tildd[p->nspin], de_tildgd[p->nspin*3];

  double C, dCdcsi, dCdzeta;
  double dzetadr[p->nspin], dcsidr[p->nspin], dcsidgr[p->nspin*3];
  int i, is;

  { /* get the PBE stuff */
    gga_type *aux2 = (p->nspin == XC_UNPOLARIZED) ? p->gga_aux2 : p->gga_aux1;
    
    gga(p->gga_aux1, rho, grho, &e_PBE, de_PBEdd, de_PBEdgd);
    
    for(is=0; is<p->nspin; is++){
      double r1[2], gr1[2], e1, de1dd[2], de1dgd[2];
      double sfac = (p->nspin == XC_UNPOLARIZED) ? 0.5 : 1.0;
      
      /* build fully polarized density and gradient */
      r1[0] = rho[is] * sfac;
      r1[1] = 0.0;
      
      for(i=0; i<3; i++){
	gr1 _(0, i) = grho _(is, i) * sfac;
	gr1 _(1, i) = 0.0;
      }

      /* call (polarized) PBE again */
      gga(aux2, r1, gr1, &e1, de1dd, de1dgd);

      if(e1 > e_PBE){
	e_til   [is] = e1;
	de_tildd[is] = de1dd[0];
	for(i=0; i<3; i++) de_tildgd _(is, i) = de1dgd _(0, i);
      }else{
	e_til   [is] = e_PBE;
	de_tildd[is] = de_PBEdd[is];
	for(i=0; i<3; i++) de_tildgd _(is, i) = de_PBEdgd _(is, i);
      }
    }
  } /* end PBE stuff */


  if(p->nspin == XC_UNPOLARIZED){
    C       = 0.53;

  }else{ /* get C(csi, zeta) */
    double gzeta[3], gzeta2, csi, a;
    
    for(i=0, gzeta2=0.0; i<3; i++){
      gzeta[i] = 2.0*(grho[0]*rho[1] - grho[1]*rho[0])/(dens*dens);
      gzeta2  += gzeta[i]*gzeta[i];
    }
    gzeta2 = max(gzeta2, MIN_GRAD*MIN_GRAD);
    
    a = 2.0*pow(3.0*M_PI*M_PI*dens, 1.0/3.0);
    csi = sqrt(gzeta2)/a;
  
    c_tpss_C(csi, zeta, &C, &dCdcsi, &dCdzeta);
    
    dzetadr[0] =  (1.0 - zeta)/dens;
    dzetadr[1] = -(1.0 + zeta)/dens;

    dcsidr [0] = -1.0/(3.0*csi); 
    dcsidr [1] = dcsidr [0];
    for(i=0; i<3; i++){
      double aux = csi*(rho[1]*grho _(0, i) - rho[0]*grho _(1, i));
      dcsidr[0] -= grho _(1, i)*aux;
      dcsidr[1] += grho _(0, i)*aux;
      
      dcsidgr _(0, i) =  rho[1]*aux;
      dcsidgr _(1, i) = -rho[0]*aux;
    }    
  }

  { /* end */
    double z2 = z*z, aux1 = 1.0 + C*z2;

    *e_PKZB    = e_PBE*aux1;
    *de_PKZBdz = dens*e_PBE*C*2.0*z;

    for(is=0; is<p->nspin; is++){
      double dCdd;

      *e_PKZB += (1.0 + C) * z2 * e_til[is] * rho[is]/dens;

      dCdd = dCdzeta*dzetadr[is] + dCdcsi*dcsidr[is];

      de_PKZBdd[is] = de_PBEdd[is]*aux1 + dCdd*dens*e_PBE*z2;
      de_PKZBdd[is]+= z2*(dCdd*e_til[is]*rho[is] +
			  (1.0 + C)*de_tildd[is]*rho[is]/dens +
			  (1.0 + C)*e_til[is]*(dens - rho[is])/dens);

      for(i=0; i<3; i++){
	double dCdgd =  dCdcsi*dcsidgr _(is, i);
	
	de_PKZBdgd _(is, i) = de_PBEdgd _(is, i)*aux1 + dCdgd*dens*e_PBE*z2;
	de_PKZBdgd _(is, i)+= z2*(dCdgd*e_til[is]*rho[is] +
				  (1.0 + C)*de_tildgd _(is, i)*rho[is]/dens);

      }

      *de_PKZBdz += (1.0 + C) * 2.0*z * e_til[is] * rho[is];
    }
  }

}


void 
mgga_c_tpss(mgga_type *p, double *rho, double *grho, double *tau,
	    double *energy, double *dedd, double *dedgd, double *dedtau)
{
  double dens, zeta;
  double gr[3], gdms, taut, tauw, z;
  double e_PKZB, de_PKZBdd[p->nspin], de_PKZBdgd[p->nspin*3], de_PKZBdz;
  int i, is;

  /* change variables */
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  /* sum tau and grho over spin */
  for(i=0; i<3; i++) gr[i] = grho _(0, i);
  taut = tau[0];

  for(is=1; is<p->nspin; is++){
    for(i=0; i<3; i++) gr[i] += grho _(is, i);
    taut += tau[is];
  }
  taut = max(taut, MIN_TAU);

  gdms = gr[0]*gr[0] + gr[1]*gr[1] + gr[2]*gr[2];
  gdms = max(MIN_GRAD*MIN_GRAD, gdms);
  
  tauw = gdms/(8.0*dens);
  z    = tauw/taut;

  /* Equation (12) */
  c_tpss_12(p, rho, grho, z, dens, zeta,
	    &e_PKZB, de_PKZBdd, de_PKZBdgd, &de_PKZBdz);

  /* Equation (11) */
  {
    double z2 = z*z, z3 = z2*z;
    double dedz;
    double dzdr, dzdgr[3], dzdtau;

    dzdr   = -z/dens;
    dzdtau = -z/taut;
    for(i=0; i<3; i++) dzdgr[i] = gr[i]/(4.0*dens*taut);

    dedz = de_PKZBdz*(1.0 + d*e_PKZB*z3) +
      d*3.0*z2*dens*e_PKZB + z3*de_PKZBdz;

    *energy = e_PKZB*(1.0 + d*e_PKZB*z3);
    for(is=0; is<p->nspin; i++){
      dedd[is]   = de_PKZBdd[is] * (1.0 + 2.0*d*e_PKZB*z3);
      dedd[is]  -= e_PKZB * d * e_PKZB * z3;
      dedd[is]  += dedz*dzdr;
    
      for(i=0; i<3; i++){
	dedgd _(is, i) = de_PKZBdgd _(is, i) * (1.0 + 2.0*d*e_PKZB*z3);
	dedgd _(is, i)+= dedz*dzdgr[i];
      }

      dedtau[is] = dedz*dzdtau; 
    }
  }

}
