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

#define XC_LDA_XC_TETER93     20   /* Teter 93 parametrization                */

static double teter_a [4] = {0.4581652932831429, 2.217058676663745,  0.7405551735357053, 0.01968227878617998 };
static double teter_ap[4] = {0.119086804055547,  0.6157402568883345, 0.1574201515892867, 0.003532336663397157};
static double teter_b [4] = {1.0000000000000000, 4.504130959426697,  1.110667363742916,  0.02359291751427506 };
static double teter_bp[4] = {0.000000000000000,  0.2673612973836267, 0.2052004607777787, 0.004200005045691381};
  

/* the functional */
void lda_c_teter93(const void *p_, const double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  double dens, zeta;
  double rs[5], aa[4], bb[4];
  double fzeta, Dfzeta;
  double Dec_Drs, Dec_Dz, ec0;

  double num, denom;
  double DnumDrs, DdenomDrs, DnumDzeta, DdenomDzeta;
  int ii;

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  /* Wigner radius */
  rs[0] = 1.0;
  rs[1] = RS(dens);
  rs[2] = rs[1]*rs[1];
  rs[3] = rs[1]*rs[2];
  rs[4] = rs[1]*rs[3];
  
  fzeta = FZETA(zeta);
  for(ii=0; ii<4; ii++){
    aa[ii] = teter_a[ii] + teter_ap[ii]*fzeta;
    bb[ii] = teter_b[ii] + teter_bp[ii]*fzeta;
  }

  /* ec */
  {
    num   = aa[0]*rs[0] + aa[1]*rs[1] + aa[2]*rs[2] + aa[3]*rs[3];
    denom = bb[0]*rs[1] + bb[1]*rs[2] + bb[2]*rs[3] + bb[3]*rs[4];
    ec0 = -num/(denom);

    if(ec != NULL)
      *ec = ec0;
  }

  /* vc */
  {
    Dfzeta      = DFZETA(zeta);

    DnumDrs     = aa[1] + 2*aa[2]*rs[1] + 3*aa[3]*rs[2];
    DdenomDrs   = bb[0] + 2*bb[1]*rs[1] + 3*bb[2]*rs[2] + 4*bb[3]*rs[3];
    Dec_Drs     = -(DnumDrs*denom - DdenomDrs*num)/(denom*denom);

    DnumDzeta   = (teter_ap[0]*rs[0] + teter_ap[1]*rs[1] + teter_ap[2]*rs[2] + teter_ap[3]*rs[3])*Dfzeta;
    DdenomDzeta = (teter_bp[0]*rs[1] + teter_bp[1]*rs[2] + teter_bp[2]*rs[3] + teter_bp[3]*rs[4])*Dfzeta;
    Dec_Dz      = -(DnumDzeta*denom - DdenomDzeta*num)/(denom*denom);

    if(vc != NULL){
      if(p->nspin == XC_UNPOLARIZED){
	vc[0] = ec0 - (rs[1]/3.0)*Dec_Drs;
      }else{
	vc[0] = ec0 - (rs[1]/3.0)*Dec_Drs - (zeta - 1.0)*Dec_Dz;
	vc[1] = ec0 - (rs[1]/3.0)*Dec_Drs - (zeta + 1.0)*Dec_Dz;
      }
    }
  }


}


const xc_func_info_type func_info_lda_xc_teter93 = {
  XC_LDA_XC_TETER93,
  XC_EXCHANGE_CORRELATION,
  "Teter 93",
  XC_FAMILY_LDA,
  "S Goedecker, M Teter, J Hutter, PRB 54, 1703 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,
  NULL,
  lda_c_teter93
};
