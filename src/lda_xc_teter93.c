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

static FLOAT teter_a [4] = {0.4581652932831429, 2.217058676663745,  0.7405551735357053, 0.01968227878617998 };
static FLOAT teter_ap[4] = {0.119086804055547,  0.6157402568883345, 0.1574201515892867, 0.003532336663397157};
static FLOAT teter_b [4] = {1.0000000000000000, 4.504130959426697,  1.110667363742916,  0.02359291751427506 };
static FLOAT teter_bp[4] = {0.000000000000000,  0.2673612973836267, 0.2052004607777787, 0.004200005045691381};
  

/* the functional */
static inline void 
func(const XC(lda_type) *p, int order, FLOAT *rs, FLOAT zeta, 
     FLOAT *zk, FLOAT *dedrs, FLOAT *dedz, 
     FLOAT *d2edrs2, FLOAT *d2edrsz, FLOAT *d2edz2)
{
  FLOAT mrs[5], aa[4], bb[4];
  FLOAT fz, dfz, d2fz;

  FLOAT num, denom, denom2, denom3;
  FLOAT DnumDrs, DdenomDrs, DnumDz, DdenomDz;
  FLOAT D2numDrs2, D2numDz2, D2numDrsz, D2denomDrs2, D2denomDz2, D2denomDrsz;
  int ii;

  /* Wigner radius */
  mrs[0] = 1.0;
  mrs[1] = rs[1];
  mrs[2] = rs[2];
  mrs[3] = mrs[1]*mrs[2];
  mrs[4] = mrs[1]*mrs[3];
  
  fz = FZETA(zeta);
  for(ii=0; ii<4; ii++){
    aa[ii] = teter_a[ii] + teter_ap[ii]*fz;
    bb[ii] = teter_b[ii] + teter_bp[ii]*fz;
  }

  num   = aa[0]*mrs[0] + aa[1]*mrs[1] + aa[2]*mrs[2] + aa[3]*mrs[3];
  denom = bb[0]*mrs[1] + bb[1]*mrs[2] + bb[2]*mrs[3] + bb[3]*mrs[4];
  *zk = -num/(denom);

  if(order < 1) return; /* nothing else to do */

  dfz       = DFZETA(zeta);
  DnumDrs   = aa[1] + 2*aa[2]*mrs[1] + 3*aa[3]*mrs[2];
  DdenomDrs = bb[0] + 2*bb[1]*mrs[1] + 3*bb[2]*mrs[2] + 4*bb[3]*mrs[3];

  DnumDz    = (teter_ap[0]*mrs[0] + teter_ap[1]*mrs[1] + teter_ap[2]*mrs[2] + teter_ap[3]*mrs[3])*dfz;
  DdenomDz  = (teter_bp[0]*mrs[1] + teter_bp[1]*mrs[2] + teter_bp[2]*mrs[3] + teter_bp[3]*mrs[4])*dfz;

  denom2 = denom*denom;

  *dedrs = -(DnumDrs*denom - DdenomDrs*num)/denom2;
  *dedz  = -(DnumDz*denom  - DdenomDz*num) /denom2;

  if(order < 2) return; /* nothing else to do */

  d2fz = D2FZETA(zeta);

  D2numDrs2   = 2*aa[2] + 3*2*aa[3]*mrs[1];
  D2denomDrs2 = 2*bb[1] + 3*2*bb[2]*mrs[1] + 4*3*bb[3]*mrs[2];

  D2numDrsz   = (teter_ap[1] + 2*teter_ap[2]*mrs[1] + 3*teter_ap[3]*mrs[2])*dfz;
  D2denomDrsz = (teter_bp[0] + 2*teter_bp[1]*mrs[1] + 3*teter_bp[2]*mrs[2] + 4*teter_bp[3]*mrs[3])*dfz;

  D2numDz2    = (teter_ap[0]*mrs[0] + teter_ap[1]*mrs[1] + teter_ap[2]*mrs[2] + teter_ap[3]*mrs[3])*d2fz;
  D2denomDz2  = (teter_bp[0]*mrs[1] + teter_bp[1]*mrs[2] + teter_bp[2]*mrs[3] + teter_bp[3]*mrs[4])*d2fz;

  denom3      = denom*denom2;

  *d2edrs2    = -((D2numDrs2*denom - D2denomDrs2*num)*denom -
		  2*DdenomDrs*(DnumDrs*denom - DdenomDrs*num))/denom3;
  *d2edz2     = -((D2numDz2*denom  - D2denomDz2*num)*denom -
		  2*DdenomDz* (DnumDz*denom  - DdenomDz*num)) /denom3;
  *d2edrsz    = -((D2numDrsz*denom + DnumDrs*DdenomDz - D2denomDrsz*num - DdenomDrs*DnumDz)*denom -
		  2*DdenomDz* (DnumDrs*denom - DdenomDrs*num))/denom3;

}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_xc_teter93) = {
  XC_LDA_XC_TETER93,
  XC_EXCHANGE_CORRELATION,
  "Teter 93",
  XC_FAMILY_LDA,
  "S Goedecker, M Teter, J Hutter, PRB 54, 1703 (1996)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,     /* init */
  NULL,     /* end  */
  work_lda, /* lda  */
};
