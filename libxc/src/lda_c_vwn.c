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

typedef struct{
  int spin_interpolation; /* 0: VWN; 1: HL */
} lda_c_vwn_params;

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
static void
init_vwn_constants(vwn_consts_type *X)
{
  int i;

  X->A[2] = -1.0/(6.0*M_PI*M_PI);
  for(i=0; i<3; i++){
    X->Q[i] = sqrt(4.0*X->c[i] - X->b[i]*X->b[i]);
  }
  X->fpp = 4.0/(9.0*(POW(2.0, 1.0/3.0) - 1));
}

static void
lda_c_vwn_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;
  lda_c_vwn_params *params;
  int func;

  assert(p->params == NULL);

  p->params = malloc(sizeof(lda_c_vwn_params));
  params = (lda_c_vwn_params *) (p->params);

  params->spin_interpolation = 0;

  func = p->info->number - XC_LDA_C_VWN;
  assert(func==0 || func==1);

  init_vwn_constants(&vwn_consts[func]);
}


void XC(lda_c_vwn_set_params)(const XC(lda_type) *p, int spin_interpolation)
{
  lda_c_vwn_params *params;

  assert(p->params != NULL);
  params = (lda_c_vwn_params *) (p->params);

  params->spin_interpolation = spin_interpolation;
}


/* Eq. (4.4) of [1] */
static void
ec_i(vwn_consts_type *X, int order, int i, FLOAT x, FLOAT *zk, FLOAT *dedrs, FLOAT *d2edrs2)
{
  FLOAT f1, f2, f3, fx, qx, xx0, drs, t1, t2, t3;
  
  /* constants */
  f1  = 2.0*X->b[i]/X->Q[i];
  f2  = X->b[i]*X->x0[i]/(X->x0[i]*X->x0[i] + X->b[i]*X->x0[i] + X->c[i]);
  f3  = 2.0*(2.0*X->x0[i] + X->b[i])/X->Q[i];

  /* a couple of handy functions */
  fx  = x*x + X->b[i]*x + X->c[i];  /* X(x) */
  qx  = atan(X->Q[i]/(2.0*x + X->b[i]));
  xx0 = x - X->x0[i];
  
  *zk = X->A[i]*(log(x*x/fx) + (f1 - f2*f3)*qx - f2*log(xx0*xx0/fx));
  
  if(order < 1) return; /* nothing else to do */

  t1 = 2.0*x + X->b[i];
  t2 = 2.0*X->c[i] + X->b[i]*x;
  t3 = t1*t1 + X->Q[i]*X->Q[i];

  drs  = X->A[i];
  drs *= -2.0*f2/xx0 + (f2*t1 + t2/x)/fx 
    - 2.0*X->Q[i]*(f1 - f2*f3)/t3;

  *dedrs = drs/(2.0*x); /* change of sqrt(rs) -> rs */

  if(order < 2) return; /* nothing else to do */
  
  *d2edrs2  = X->A[i];
  *d2edrs2 *= -f2*t1*t1/(fx*fx) - t1*t2/(x*fx*fx) + 2.0*f2/fx
    + X->b[i]/(x*fx) - t2/(x*x*fx) + 8.0*(f1 - f2*f3)*X->Q[i]*t1/(t3*t3)
    + 2.0*f2/(xx0*xx0);

  *d2edrs2  = ((*d2edrs2)*x - drs)/(2.0*x*x);
  *d2edrs2 /= 2.0*x;
}

/* the functional */
static inline void 
func(const XC(lda_type) *p, int order, FLOAT *rs, FLOAT zeta, 
     FLOAT *zk, FLOAT *dedrs, FLOAT *dedz, 
     FLOAT *d2edrs2, FLOAT *d2edrsz, FLOAT *d2edz2)
{
  int func;
  vwn_consts_type *X;
  lda_c_vwn_params *params;

  func = p->info->number - XC_LDA_C_VWN;
  assert(func==0 || func==1);

  assert(p->params != NULL);
  params = (lda_c_vwn_params *) (p->params);

  X = &vwn_consts[func];

  ec_i(X, order, 0, rs[0], zk, dedrs, d2edrs2);
  
  if(p->nspin==XC_POLARIZED){
    FLOAT ec1, ec2, ec3, vc1, vc2, vc3, fc1, fc2, fc3;
    FLOAT z3, z4, t1, dt1, d2t1, t2, dt2, d2t2, fz, dfz, d2fz;
    
    /* store paramagnetic values */
    ec1 = *zk;
    if(order >= 1) vc1 = *dedrs;
    if(order >= 2) fc1 = *d2edrs2;
    
    ec_i(X, order, 1, rs[0], &ec2, &vc2, &fc2);
    ec_i(X, order, 2, rs[0], &ec3, &vc3, &fc3);
    
    fz  = FZETA(zeta);
    if(params->spin_interpolation == 1){
      t1 = 0.0;
      t2 = fz;
    }else{
      z3  = POW(zeta, 3);
      z4  = z3*zeta;
      t1  = (fz/X->fpp)*(1.0 - z4);
      t2  = fz*z4;
    }

    *zk =  ec1 +  ec3*t1 + (ec2 -  ec1)*t2;

    if(order < 1) return; /* nothing else to do */

    dfz  = DFZETA(zeta);
    if(params->spin_interpolation == 1){
      dt1 = 0.0;
      dt2 = dfz;
    }else{
      dt1  = dfz*(1.0 - z4) - 4.0*fz*z3;
      dt1 /= X->fpp;
      dt2  = dfz*z4 + 4.0*fz*z3;
    }

    *dedrs = vc1 + vc3* t1 + (vc2 - vc1)* t2;
    *dedz  =       ec3*dt1 + (ec2 - ec1)*dt2;

    if(order < 2) return; /* nothing else to do */

    d2fz  = D2FZETA(zeta);
    if(params->spin_interpolation == 1){
      dt1 = 0.0;
      dt2 = d2fz;
    }else{
      d2t1  = d2fz*(1.0 - z4) - 8.0*dfz*z3 - 4.0*3.0*fz*zeta*zeta;
      d2t1 /= X->fpp;
      d2t2  = d2fz*z4 + 8.0*dfz*z3 + 4.0*3.0*fz*zeta*zeta;
    }

    *d2edrs2 = fc1 + fc3*  t1 + (fc2 - fc1)*  t2;
    *d2edrsz =       vc3* dt1 + (vc2 - vc1)* dt2;
    *d2edz2  =       ec3*d2t1 + (ec2 - ec1)*d2t2;
  }
	
}

#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_c_vwn) = {
  XC_LDA_C_VWN,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair",
  XC_FAMILY_LDA,
  "SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  lda_c_vwn_init,
  NULL,
  work_lda
};

const XC(func_info_type) XC(func_info_lda_c_vwn_rpa) = {
  XC_LDA_C_VWN_RPA,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair (parametrization of the RPA energy)",
  XC_FAMILY_LDA,
  "SH Vosko, L Wilk, and M Nusair, Can. J. Phys. 58, 1200 (1980)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  lda_c_vwn_init,
  NULL,
  work_lda 
};

