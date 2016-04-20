/*
 Copyright (C) 2014 Orbital-free DFT group at University of Florida, USA

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/***********************************************************************
  Exchange and correlation free energy density and potential as 
  parametrized by 
    Valentin V. Karasiev, Travis Sjostrom, James Dufty, and S. B. Trickey
  Ported to C and libxc by Lazaro Calderin and Miguel Marques
************************************************************************/

#define XC_LDA_XC_KSDT    259    /* Karasiev et al. parametrization */

static const FLOAT ksdt_a[6] =
  {0.750, 3.043630, -0.0922700, 1.703500, 8.310510, 5.11050};
static const FLOAT ksdt_b[2][5] = {               /* b5 = Sqrt[3/2]/(lambda)*b3 */
  {0.2839970,  48.9321540, 0.3709190, 61.0953570, 0.871837422702767684673873513724},
  {0.3290010, 111.5983080, 0.5370530,105.0866630, 1.26233194679913807935662124247}
};
static const FLOAT ksdt_c[2][3] = {
  {0.8700890, 0.1930770, 2.4146440},
  {0.8489300, 0.1679520, 0.0888200}
};
static const FLOAT ksdt_d[2][5] = {
  {0.5798240,  94.5374540,  97.8396030,  59.9399990, 24.3880370},
  {0.5513300, 180.2131590, 134.4862310, 103.8616950, 17.7507100}
};
static const FLOAT ksdt_e[2][5] = {
  {0.2120360, 16.7312490, 28.4857920,  34.0288760, 17.2355150},
  {0.1531240, 19.5439450, 43.4003370, 120.2551450, 15.6628360}
};

typedef struct{
  FLOAT T;
} lda_xc_ksdt_params;


static void 
lda_xc_ksdt_init(XC(func_type) *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_xc_ksdt_params));
  XC(lda_xc_ksdt_set_params)(p, 0.0);
}

void 
XC(lda_xc_ksdt_set_params)(XC(func_type) *p, FLOAT T)
{
  lda_xc_ksdt_params *params;

  assert(p != NULL && p->params != NULL);
  params = (lda_xc_ksdt_params *) (p->params);

  params->T  = T;
}

void
ksdt_fxc(int ispin, int order, FLOAT t, FLOAT *rs, FLOAT *fxc, FLOAT *dfxcdt, 
	 FLOAT *dfxcdrs, FLOAT *d2fxdt2, FLOAT *d2fxcrsdt, FLOAT *d2fxcdrs2)
{
  const FLOAT lambda = 0.521061761197848019684674268560; /* pow(4/(9*pi), 1.0/3.0) */
  const FLOAT a0     = 0.610887057710857191300739313402; /* 1/(pi*lambda)          */

  FLOAT omega;
  FLOAT t2, t3, t4, sqrtt, tanht, tanhsqrt;
  FLOAT a_num, a_den, aa, b_num, b_den, bb, d_num, d_den, dd, e_num, e_den, ee, c_num, cc;
  FLOAT f1, fxc_num, fxc_den;

  FLOAT dtanht, dtanhsqrt;
  FLOAT da_num, da_den, daa, db_num, db_den, dbb, dd_num, dd_den, ddd, de_num, de_den, dee, dc_num, dcc;
  FLOAT df1, dfxc_numdt, dfxc_dendt, dfxc_numdrs, dfxc_dendrs;


  omega = (ispin == 0) ? 1.0 : M_CBRT2;

  t2 = t*t; t3 = t2*t; t4 = t3*t;

  sqrtt = SQRT(t);
  tanht = tanh(1.0/t);
  tanhsqrt = tanh(1.0/sqrtt);

  a_num = ksdt_a[0] + ksdt_a[1]*t2 + ksdt_a[2]*t3 + ksdt_a[3]*t4;
  a_den = 1.0 + ksdt_a[4]*t2 + ksdt_a[5]*t4;
  aa    = a0*tanht*a_num/a_den;
  
  b_num = ksdt_b[ispin][0] + ksdt_b[ispin][1]*t2 + ksdt_b[ispin][2]*t4;
  b_den = 1.0 + ksdt_b[ispin][3]*t2 + omega*ksdt_b[ispin][4]*t4;
  bb    = tanhsqrt*b_num/b_den;

  d_num = ksdt_d[ispin][0] + ksdt_d[ispin][1]*t2 + ksdt_d[ispin][2]*t4;
  d_den = 1.0 + ksdt_d[ispin][3]*t2 + ksdt_d[ispin][4]*t4;
  dd    = tanhsqrt*d_num/d_den;

  e_num = ksdt_e[ispin][0] + ksdt_e[ispin][1]*t2 + ksdt_e[ispin][2]*t4;
  e_den = 1.0 + ksdt_e[ispin][3]*t2 + ksdt_e[ispin][4]*t4;
  ee    = tanht*e_num/e_den;

  c_num = ksdt_c[ispin][0] + ksdt_c[ispin][1]*EXP(-ksdt_c[ispin][2]/t);
  cc    = c_num*ee;

  f1      = -1.0/rs[1];
  fxc_num = omega*aa + bb*rs[0] + cc*rs[1];
  fxc_den = 1.0 + dd*rs[0] + ee*rs[1];
  *fxc    = f1*fxc_num/fxc_den;

  if(order < 1) return;

  /* derivatives */
  dtanht    = (tanht*tanht - 1.0)/t2;
  dtanhsqrt = (tanhsqrt*tanhsqrt - 1.0)/(2.0*t*sqrtt);

  da_num = ksdt_a[1]*2.0*t + ksdt_a[2]*3.0*t2 + ksdt_a[3]*4.0*t3;
  da_den = ksdt_a[4]*2.0*t + ksdt_a[5]*4.0*t3;
  daa    = a0*DFRACTION(tanht*a_num, dtanht*a_num + tanht*da_num, a_den, da_den);

  db_num = ksdt_b[ispin][1]*2.0*t + ksdt_b[ispin][2]*4.0*t3;
  db_den = ksdt_b[ispin][3]*2.0*t + omega*ksdt_b[ispin][4]*4.0*t3;
  dbb    = DFRACTION(tanhsqrt*b_num, dtanhsqrt*b_num + tanhsqrt*db_num, b_den, db_den);

  dd_num = ksdt_d[ispin][1]*2.0*t + ksdt_d[ispin][2]*4.0*t3;
  dd_den = ksdt_d[ispin][3]*2.0*t + ksdt_d[ispin][4]*4.0*t3;
  ddd    = DFRACTION(tanhsqrt*d_num, dtanhsqrt*d_num + tanhsqrt*dd_num, d_den, dd_den);

  de_num = ksdt_e[ispin][1]*2.0*t + ksdt_e[ispin][2]*4.0*t3;
  de_den = ksdt_e[ispin][3]*2.0*t + ksdt_e[ispin][4]*4.0*t3;
  dee    = DFRACTION(tanht*e_num, dtanht*e_num + tanht*de_num, e_den, de_den);

  dc_num = ksdt_c[ispin][1]*ksdt_c[ispin][2]*exp(-ksdt_c[ispin][2]/t)/t2;
  dcc     = dc_num*ee + c_num*dee;

  dfxc_numdt = omega*daa + dbb*rs[0] + dcc*rs[1];
  dfxc_dendt = ddd*rs[0] + dee*rs[1];
  *dfxcdt    = DFRACTION(f1*fxc_num, f1*dfxc_numdt, fxc_den, dfxc_dendt);

  df1 = -f1/rs[1];
  dfxc_numdrs = bb/(2.0*rs[0]) + cc;
  dfxc_dendrs = dd/(2.0*rs[0]) + ee;

  *dfxcdrs    = DFRACTION(f1*fxc_num, df1*fxc_num + f1*dfxc_numdrs, fxc_den, dfxc_dendrs);

  if(order < 2) return;

}

void
ksdt_alpha(int order, FLOAT t, FLOAT *rs, FLOAT *alpha, FLOAT *dalphadt, FLOAT *dalphadrs)
{
  const FLOAT ksdt_g[] = {2.0/3.0, -0.0139261, 0.183208};
  const FLOAT ksdt_l[] = {1.064009, 0.572565};

  FLOAT ll, gg_num, gg_den, gg, aux;
  FLOAT dlldt, dlldrs, dggdrs, dauxdt, dauxdrs;

  ll = ksdt_l[0] + ksdt_l[1]*t*rs[0];

  gg_num = ksdt_g[0] + ksdt_g[1]*rs[1];
  gg_den = 1.0 + ksdt_g[2]*rs[1];
  gg     = gg_num/gg_den;

  aux    = EXP(-t*ll);

  *alpha = 2.0 - gg*aux;

  if(order < 1) return;

  dlldt  = ksdt_l[1]*rs[0];
  dlldrs = ksdt_l[1]*t/(2.0*rs[0]);

  dggdrs = DFRACTION(gg_num, ksdt_g[1], gg_den, ksdt_g[2]);

  dauxdt  = -(ll + t*dlldt)*aux;
  dauxdrs = -t*dlldrs*aux;

  *dalphadt  = -gg*dauxdt;
  *dalphadrs = -(dggdrs*aux + gg*dauxdrs);
}


void
ksdt_phi(int order, FLOAT zeta, FLOAT alpha, 
	 FLOAT *phi, FLOAT *dphidz, FLOAT *dphidalpha)
{
  FLOAT fzetafactor, opz, omz, opza, omza;

  fzetafactor = POW(2.0, alpha) - 2.0;

  opz = 1.0 + zeta;
  omz = 1.0 - zeta;
  
  opza = POW(opz, alpha);
  omza = POW(omz, alpha);

  *phi = (opza + omza - 2.0)/fzetafactor;

  if(order < 1) return;

//  *dphidz = alpha*(opza/opz - omza/omz)/fzetafactor;
  *dphidz = alpha*( POW(opz,alpha-1.0) - POW(omz,alpha-1.0) )/fzetafactor;

  if( omz != 0.0 ){
  *dphidalpha = DFRACTION(opza + omza - 2.0, opza*LOG(opz) + omza*LOG(omz), 
			  fzetafactor, POW(2.0, alpha)*LOG(2.0));
  }else{
  *dphidalpha = DFRACTION(opza + omza - 2.0, opza*LOG(opz),
                          fzetafactor, POW(2.0, alpha)*LOG(2.0));
  } 
}

/* the functional */
void 
XC(lda_xc_ksdt2)(const XC(func_type) *p, XC(lda_work_t) *r)
{
  FLOAT temp, tr_factor, tt, dtdrs;
  FLOAT fxc0, dfxc0dt, dfxc0drs, d2fxc0dt2, d2fxc0drst, d2fxc0drs2;
  FLOAT fxc1, dfxc1dt, dfxc1drs, d2fxc1dt2, d2fxc1drst, d2fxc1drs2;
  FLOAT alpha, dalphadt, dalphadrs;
  FLOAT phi, dphidz, dphidalpha;

  assert(p->params != NULL);
  temp = max(((lda_xc_ksdt_params *) (p->params))->T, 1e-8);

  //tr_factor = 2.0/CBRT(3.0*M_PI*M_PI);
  //tr_factor *= tr_factor*temp/(RS_FACTOR*RS_FACTOR);

  tr_factor  = CBRT(4.0/(9.0*M_PI));
  tr_factor *= tr_factor*2.0*temp; 


  tt = tr_factor*r->rs[2];

  ksdt_fxc(0, r->order, tt, r->rs, 
	   &fxc0, &dfxc0dt, &dfxc0drs, &d2fxc0dt2, &d2fxc0drst, &d2fxc0drs2);

  if(p->nspin == XC_UNPOLARIZED){
    r->zk = fxc0;
  }else{
    ksdt_fxc(1, r->order, tt/(M_CBRT2*M_CBRT2), r->rs, 
	     &fxc1, &dfxc1dt, &dfxc1drs, &d2fxc1dt2, &d2fxc1drst, &d2fxc1drs2);

    ksdt_alpha(r->order, tt, r->rs, &alpha, &dalphadt, &dalphadrs);
    ksdt_phi(r->order, r->zeta, alpha, &phi, &dphidz, &dphidalpha);

    r->zk = fxc0 + (fxc1 - fxc0)*phi;
  }

  if(r->order < 1) return;

  dtdrs = tr_factor*2.0*r->rs[1];

  if(p->nspin == XC_UNPOLARIZED){
    r->dedrs = dfxc0drs + dfxc0dt*dtdrs;
  }else{
    dfxc1dt /= M_CBRT2*M_CBRT2;

    r->dedrs = dfxc0drs + (dfxc1drs - dfxc0drs)*phi +
      (dfxc0dt + (dfxc1dt - dfxc0dt)*phi)*dtdrs +
      (fxc1 - fxc0)*dphidalpha*(dalphadrs + dalphadt*dtdrs);
    r->dedz  = (fxc1 - fxc0)*dphidz;
  }
}

#define func XC(lda_xc_ksdt2)
#include "work_lda.c"

const XC(func_info_type) XC(func_info_lda_xc_ksdt) = {
  XC_LDA_XC_KSDT,
  XC_EXCHANGE_CORRELATION,
  "Karasiev, Sjostrom, Dufty & Trickey",
  XC_FAMILY_LDA,
  {&xc_ref_Karasiev2014_076403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 0.0, 0.0, 1e-32,
  lda_xc_ksdt_init,     /* init */
  NULL,     /* end  */
  work_lda,  /* lda  */
  NULL,
  NULL
};
