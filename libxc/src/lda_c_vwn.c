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

#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 LDA parametrization of Vosko, Wilk & Nusair
************************************************************************/

#define XC_LDA_C_VWN      7   /* Vosko, Wilk, & Nussair       */
#define XC_LDA_C_VWN_RPA  8   /* Vosko, Wilk, & Nussair (RPA) */

/* some constants         e_c^P      e_c^F      alpha_c */
typedef struct {
  FLOAT  A[3]; /* e_c^P, e_c^F, alpha_c */
  FLOAT  b[3];
  FLOAT  c[3];
  FLOAT x0[3];
  FLOAT  Q[3];
  FLOAT  fpp;
} vwn_consts_type;

/* These numbers are taken from the original reference, but divided by
     two to convert from Rydbergs to Hartrees */
static vwn_consts_type vwn_consts[2] = {
  /* VWN parametrization of the correlation energy */
  {
    { 0.0310907, 0.01554535,  0.0      }, /*  A */
    { 3.72744,   7.06042,     1.13107  }, /*  b */
    {12.9352,   18.0578,     13.0045   }, /*  c */
    {-0.10498,  -0.32500,    -0.0047584}, /* x0 */
    { 0.0,       0.0,         0.0      }, /*  Q */
    0.0 /* fpp */
  },
  /* VWN RPA */
  {
    { 0.0310907, 0.01554535,  0.0      }, /*  A */
    {13.0720,   20.1231,      1.06835  }, /*  b */
    {42.7198,  101.578,      11.4813   }, /*  c */
    {-0.409286, -0.743294,   -0.228344 }, /* x0 */
    { 0.0,       0.0,         0.0      }, /*  Q */
    0.0 /* fpp */
  }
};

/* initialization */
void init_vwn_constants(vwn_consts_type *X)
{
  int i;

  X->A[2] = -1.0/(6.0*M_PI*M_PI);
  for(i=0; i<3; i++){
    X->Q[i] = sqrt(4.0*X->c[i] - X->b[i]*X->b[i]);
  }
  X->fpp = 4.0/(9.0*(POW(2.0, 1.0/3.0) - 1));
}

static void lda_c_vwn_init(void *p_)
{
  init_vwn_constants(&vwn_consts[0]);
}

static void lda_c_vwn_rpa_init(void *p_)
{
  init_vwn_constants(&vwn_consts[1]);
}

/* Eq. (4.4) of [1] */
void ec_i(vwn_consts_type *X, int i, FLOAT x, FLOAT *ec, FLOAT *decdrs)
{
  FLOAT f1, f2, f3, fx, qx, xx0, tx, tt;
  
  f1  = 2.0*X->b[i]/X->Q[i];
  f2  = X->b[i]*X->x0[i]/(X->x0[i]*X->x0[i] + X->b[i]*X->x0[i] + X->c[i]);
  f3  = 2.0*(2.0*X->x0[i] + X->b[i])/X->Q[i];
  fx  = x*x + X->b[i]*x + X->c[i];  /* X(x) */
  qx  = atan(X->Q[i]/(2.0*x + X->b[i]));
  xx0 = x - X->x0[i];
  
  *ec = X->A[i]*(log(x*x/fx) + f1*qx - f2*(log(xx0*xx0/fx) + f3*qx));
  
  tx  = 2.0*x + X->b[i];
  tt  = tx*tx + X->Q[i]*X->Q[i];
  *decdrs = X->A[i]*(2.0/x - tx/fx - 4.0*X->b[i]/tt -
		     f2*(2.0/xx0 - tx/fx - 4.0*(2.0*X->x0[i] + X->b[i])/tt));
}

/* the functional */
void lda_c_vwn(const void *p_, const FLOAT *rho, FLOAT *ec, FLOAT *vc, FLOAT *fc)
{
  const xc_lda_type *p = (xc_lda_type *)p_;

  FLOAT dens, zeta;
  FLOAT rs[2], ec_1, dec_1;
  int func;
  vwn_consts_type *X;

  assert(p!=NULL);

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  func = p->info->number - XC_LDA_C_VWN;
  assert(func==0 || func==1);
  X = &vwn_consts[func];

  /* Wigner radius */
  rs[1] = RS(dens);     /* rs          */
  rs[0] = sqrt(rs[1]);  /* sqrt(rs)    */

  ec_i(X, 0, rs[0], &ec_1, &dec_1);
  
  if(p->nspin == XC_UNPOLARIZED){
    *ec   = ec_1;
    vc[0] = ec_1 - rs[0]/6.0*dec_1;
  }else{
    FLOAT ec_2, ec_3, dec_2, dec_3, fz, dfz, decdx, decdz;
    FLOAT t1, t2, z3, z4;
    
    ec_i(X, 1, rs[0], &ec_2, &dec_2);
    ec_i(X, 2, rs[0], &ec_3, &dec_3);
    
    fz  =  FZETA(zeta);
    dfz = DFZETA(zeta);
    
    z3 = POW(zeta, 3);
    z4 = z3*zeta;
    t1 = (fz/X->fpp)*(1.0 - z4);
    t2 = fz*z4;
    
    *ec   =  ec_1 +  ec_3*t1 + ( ec_2 -  ec_1)*t2;
    decdx = dec_1 + dec_3*t1 + (dec_2 - dec_1)*t2;
    decdz = (ec_3/X->fpp)*(dfz*(1.0 - z4) - 4.0*fz*z3) +
      (ec_2 - ec_1)*(dfz*z4 + 4.0*fz*z3);
    
    t1 = *ec - rs[0]/6.0*decdx;
    
    vc[0] = t1 + (1.0 - zeta)*decdz;
    vc[1] = t1 - (1.0 + zeta)*decdz;
  }
	
}

const xc_func_info_type func_info_lda_c_vwn = {
  XC_LDA_C_VWN,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair",
  XC_FAMILY_LDA,
  "SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_c_vwn_init,
  NULL,
  lda_c_vwn
};

const xc_func_info_type func_info_lda_c_vwn_rpa = {
  XC_LDA_C_VWN_RPA,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair (parametrization of the RPA energy)",
  XC_FAMILY_LDA,
  "SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_c_vwn_rpa_init,
  NULL,
  lda_c_vwn 
};

