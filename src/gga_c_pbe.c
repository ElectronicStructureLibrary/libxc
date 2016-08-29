/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

/************************************************************************
 Implements Perdew, Burke & Ernzerhof Generalized Gradient Approximation
 correlation functional.

 I based this implementation on a routine from L.C. Balbas and J.M. Soler
************************************************************************/

#define XC_GGA_C_PBE          130 /* Perdew, Burke & Ernzerhof correlation              */
#define XC_GGA_C_PBE_SOL      133 /* Perdew, Burke & Ernzerhof correlation SOL          */
#define XC_GGA_C_XPBE         136 /* xPBE reparametrization by Xu & Goddard             */
#define XC_GGA_C_PBE_JRGX     138 /* JRGX reparametrization by Pedroza, Silva & Capelle */
#define XC_GGA_C_RGE2         143 /* Regularized PBE                                    */
#define XC_GGA_C_APBE         186 /* mu fixed from the semiclassical neutral atom       */
#define XC_GGA_C_SPBE          89 /* PBE correlation to be used with the SSB exchange   */
#define XC_GGA_C_REGTPSS       83 /* Regularized TPSS correlation (ex-VPBE)             */
#define XC_GGA_C_ZPBESOL       63 /* spin-dependent gradient correction to PBEsol       */
#define XC_GGA_C_PBEINT        62 /* PBE for hybrid interfaces                          */
#define XC_GGA_C_ZPBEINT       61 /* spin-dependent gradient correction to PBEint       */
#define XC_GGA_C_PBELOC       246 /* Semilocal dynamical correlation                    */
#define XC_GGA_C_BCGP          39 /* Burke, Cancio, Gould, and Pittalis                 */
#define XC_GGA_C_PBEFE        258 /* PBE for formation energies                         */
#define XC_GGA_C_PBE_MOL      272 /* Del Campo, Gazquez, Trickey and Vela (PBE-like)    */
#define XC_GGA_C_SG4          534 /* Semiclassical GGA at fourth order                  */

typedef struct{
  FLOAT beta;
} gga_c_pbe_params;


static void gga_c_pbe_init(XC(func_type) *p)
{
  static const FLOAT beta[]  = {
    0.06672455060314922,                /*  0: original PBE              */
    0.046,                              /*  1: PBE sol                   */
    0.089809,                           /*  2: xPBE                      */
    3.0*10.0/(81.0*M_PI*M_PI),          /*  3: PBE_JRGX                  */
    0.053,                              /*  4: RGE2                      */
    3.0*0.260/(M_PI*M_PI),              /*  5: APBE (C)                  */
    0.06672455060314922,                /*  6: sPBE                      */
    0.0,                                /*  7: vPBE this is calculated   */
    0.046,                              /*  8: zPBEsol                   */
    0.052,                              /*  9: PBEint                    */
    0.052,                              /* 10: zPBEint                   */
    0.0,                                /* 11: PBEloc this is calculated */
    0.06672455060314922,                /* 12: BCGP                      */
    0.043,                              /* 13: PBEfe                     */
    0.08384,                            /* 14: PBEmol                    */
    0.0                                 /* 15: SG14                      */
  };

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_c_pbe_params));

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PW_MOD, p->nspin);

  switch(p->info->number){
  case XC_GGA_C_PBE:      p->func =  0; break;    
  case XC_GGA_C_PBE_SOL:  p->func =  1; break;
  case XC_GGA_C_XPBE:     p->func =  2; break;
  case XC_GGA_C_PBE_JRGX: p->func =  3; break;
  case XC_GGA_C_RGE2:     p->func =  4; break;
  case XC_GGA_C_APBE:     p->func =  5; break;
  case XC_GGA_C_SPBE:     p->func =  6; break;
  case XC_GGA_C_REGTPSS:  p->func =  7; break;
  case XC_GGA_C_ZPBESOL:  p->func =  8; break;
  case XC_GGA_C_PBEINT:   p->func =  9; break;
  case XC_GGA_C_ZPBEINT:  p->func = 10; break;
  case XC_GGA_C_PBELOC:   p->func = 11; break;
  case XC_GGA_C_BCGP:     p->func = 12; break;
  case XC_GGA_C_PBEFE:    p->func = 13; break;
  case XC_GGA_C_PBE_MOL:  p->func = 14; break;
  case XC_GGA_C_SG4:      p->func = 15; break;
    /* func = 99 is reserved for the worker routine of MGGA_C_SCAN */
  default:
    fprintf(stderr, "Internal error in gga_c_pbe\n");
    exit(1);
  }

  XC(gga_c_pbe_set_params)(p, beta[p->func]);
}


void 
XC(gga_c_pbe_set_params)(XC(func_type) *p, FLOAT beta)
{
  gga_c_pbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_c_pbe_params *) (p->params);

  params->beta = beta;
}


static inline void
bcgp_pt(int order, FLOAT tt, FLOAT *tp, FLOAT *dtpdtt, FLOAT *d2tpdtt2)
{
  const FLOAT cac = 1.467;
  const FLOAT tau = 4.5;
  FLOAT num, den, P, P_2, dP, d2P;

  num = tau + tt;
  den = tau + cac*tt;
  
  P   = num/den; 
  P_2 = SQRT(P);

  *tp = tt * P_2;

  if(order < 1) return;

  dP = DFRACTION(num, 1.0, den, cac);
  *dtpdtt = P_2 + tt*dP/(2.0*P_2);

  if(order < 2) return;

  d2P = D2FRACTION(num, 1.0, 0.0, den, cac, 0.0);
  *d2tpdtt2 = (2.0*dP + tt*d2P - tt*dP*dP/(2.0*P))/(2.0*P_2);
}

static inline void 
pbe_eq8(int order, FLOAT beta, FLOAT gamm, FLOAT ecunif, FLOAT phi, 
	FLOAT *A, FLOAT *dbeta, FLOAT *dec, FLOAT *dphi,
	FLOAT *dec2, FLOAT *decphi, FLOAT *dphi2)
{
  FLOAT phi3, f1, df1dphi, d2f1dphi2, f2, f3, dx, d2x;

  phi3 = POW(phi, 3);
  f1   = ecunif/(gamm*phi3);
  f2   = EXP(-f1);
  f3   = f2 - 1.0;

  *A   = beta/(gamm*f3);

  if(order < 1) return;

  df1dphi = -3.0*f1/phi;
  dx      = (*A)*f2/f3;

  *dbeta  = 1.0/(gamm*f3);
  *dec    = dx/(gamm*phi3);
  *dphi   = dx*df1dphi;

  if(order < 2) return;

  d2f1dphi2 = -4.0*df1dphi/phi;
  d2x       = dx*(2.0*f2 - f3)/f3;
  *dphi2    = d2x*df1dphi*df1dphi + dx*d2f1dphi2;
  *decphi   = (d2x*df1dphi*f1 + dx*df1dphi)/ecunif;
  *dec2     = d2x/(gamm*gamm*phi3*phi3);
}


static inline void 
pbe_eq7(int order, int func, FLOAT beta, FLOAT gamm, FLOAT phi, FLOAT t, FLOAT A, FLOAT B,
	FLOAT *H, FLOAT *dbeta, FLOAT *dphi, FLOAT *dt, FLOAT *dA,
	FLOAT *d2phi, FLOAT *d2phit, FLOAT *d2phiA, FLOAT *d2t2, FLOAT *d2tA, FLOAT *d2A2)
{
  FLOAT alpha, t2, phi3, f1, f2, f3;
  FLOAT df2dbeta, df1dt, df2dt, df1dA, df2dA;
  FLOAT d2f1dt2, d2f2dt2, d2f2dA2, d2f1dtA, d2f2dtA;
  FLOAT ff, dffdt, dffdphi;

  t2   = t*t;
  phi3 = POW(phi, 3);

  if(func == 99){ /* XC_GGA_C_SCAN_E0 */
    f1 = 1.0 + 4.0*A*t2;
    f3 = POW(f1, -1.0/4.0);
    f2 = beta*(1.0 - f3)/(gamm*A);
  }else{
    f1 = t2 + B*A*t2*t2;
    f3 = 1.0 + A*f1;
    f2 = beta*f1/(gamm*f3);
  }

  if(func == 8 || func == 10 || func == 15){ /* zPBEsol, zPBEint, and SG4 */
    if(func == 15)
      alpha = 0.8;
    else
      alpha = (func == 8) ? 4.8 : 2.4;
    ff = POW(phi, alpha*t*t2);
  }else
    ff = 1.0;

  *H = ff*gamm*phi3*LOG(1.0 + f2);

  if(order < 1) return;

  if(func == 8 || func == 10 || func == 15){ /* zPBEsol, zPBEint, and SG4 */
    dffdphi = alpha*t*t2*ff/phi;
    dffdt   = 3.0*alpha*t2*ff*LOG(phi);
  }else{
    dffdphi = 0.0;
    dffdt   = 0.0;
  }

  *dphi  = 3.0*(*H)/phi + (*H)*dffdphi/ff;
  
  df2dbeta = f2/beta;
  *dbeta   = ff*gamm*phi3*df2dbeta/(1.0 + f2);

  if(func == 99){ /* XC_GGA_C_SCAN_E0 */
    df2dt    = 2.0*beta*t*f3/(gamm*f1);
  }else{
    df1dt    = t*(2.0 + 4.0*B*A*t2);
    df2dt    = beta/(gamm*f3*f3) * df1dt;
  }

  *dt      = ff*gamm*phi3*df2dt/(1.0 + f2) + (*H)*dffdt/ff;

  if(func == 99){ /* XC_GGA_C_SCAN_E0 */
    df2dA    = -beta/(gamm*A*A) * (1.0 - (1 + 5.0*A*t2)*f3/f1);
  }else{
    df1dA    = B*t2*t2;
    df2dA    = beta/(gamm*f3*f3) * (df1dA - f1*f1);
  }
    *dA      = ff*gamm*phi3*df2dA/(1.0 + f2);

  if(order < 2) return;

  *d2phi  = 2.0*(*dphi)/phi;
  *d2phit = 3.0*(*dt)/phi;
  *d2phiA = 3.0*(*dA)/phi;

  d2f1dt2 = 2.0 + 4.0*3.0*B*A*t2;
  d2f2dt2 = beta/(gamm*f3*f3) * (d2f1dt2 - 2.0*A/f3*df1dt*df1dt);
  *d2t2   = gamm*phi3*(d2f2dt2*(1.0 + f2) - df2dt*df2dt)/((1.0 + f2)*(1.0 + f2));

  d2f1dtA = 4.0*B*t*t2;
  d2f2dtA = beta/(gamm*f3*f3) * 
    (d2f1dtA - 2.0*df1dt*(f1 + A*df1dA)/f3);
  *d2tA   = gamm*phi3*(d2f2dtA*(1.0 + f2) - df2dt*df2dA)/((1.0 + f2)*(1.0 + f2));

  d2f2dA2 = beta/(gamm*f3*f3*f3) *(-2.0)*(2.0*f1*df1dA - f1*f1*f1 + A*df1dA*df1dA);
  *d2A2   = gamm*phi3*(d2f2dA2*(1.0 + f2) - df2dA*df2dA)/((1.0 + f2)*(1.0 + f2));
}


/* Calculates beta(r) = 0.066725 (1 + 0.1 r_s) / (1 + 0.1778 r_s) */
void
XC(beta_Hu_Langreth) (FLOAT rs, int order, FLOAT *b, FLOAT *dbdrs, FLOAT *d2bdrs2)
{
  /* in the paper we have beta_a = 0.066725 */
  const static FLOAT beta_a = 0.066724550603149220, beta_b = 0.1, beta_c = 0.1778;

  FLOAT num, den, dnum, dden;
  
  den = 1.0 + beta_c*rs;
  num = beta_a * (1.0 + beta_b*rs);
  *b =  num / den;

  if(order < 1) return;

  dden = beta_c;
  dnum = beta_a*beta_b;

  *dbdrs = DFRACTION(num, dnum, den, dden);

  if(order < 2) return;

  *d2bdrs2 = D2FRACTION(num, dnum, 0.0, den, dden, 0.0);
}


inline void 
XC(gga_c_pbe_func) (const XC(func_type) *p, XC(gga_work_c_t) *r)
{
  /* parameters for beta of PBEloc */
  static FLOAT pbeloc_b0 = 0.0375, pbeloc_a = 0.08;
  /* parameters for beta of SG4, b0 = 3 mu^MGE2/pi^2, mu^MGE2 = 0.26 */
  static FLOAT sg4_b0 = 3.0*0.262/(M_PI*M_PI), sg4_sigma = 0.07; 

  FLOAT cnst_beta, phi, tt, tp, dtpdtt, d2tpdtt2;

  FLOAT A, dAdbeta, dAdec, dAdphi, d2Adec2, d2Adecphi, d2Adphi2;
  FLOAT H, dHdbeta, dHdphi, dHdt, dHdA, d2Hdphi2, d2Hdphit, d2HdphiA, d2Hdt2, d2HdtA, d2HdA2;
  FLOAT dfdbeta, dfdphi, dfdec, dfdt, dtdrs, dtdxt, dtdphi, dphidz;
  FLOAT d2fdphi2, d2fdphit, d2fdphiec, d2fdt2, d2fdtec, d2fdec2;
  FLOAT d2tdrs2, d2tdrsxt, d2tdphi2, d2tdrsphi, d2tdxtphi, d2phidz2;
  FLOAT B;

  XC(lda_work_t) pw;
  FLOAT tconv, auxp, auxm, beta, gamm, beta_den, dbetadrs, d2betadrs2;

  assert(p->params != NULL);
  cnst_beta = ((gga_c_pbe_params *) (p->params))->beta;

  pw.order = r->order;
  pw.rs[0] = SQRT(r->rs);
  pw.rs[1] = r->rs;
  pw.rs[2] = r->rs*r->rs;
  pw.zeta  = r->zeta;

  XC(lda_c_pw_func)(p->func_aux[0], &pw);

  tconv = 4.0*M_CBRT2;

  auxp = CBRT(1.0 + r->zeta);
  auxm = CBRT(1.0 - r->zeta);

  phi  = 0.5*(auxp*auxp + auxm*auxm);
  tt   = r->xt/(tconv*phi*pw.rs[0]);

  if(p->info->number == XC_GGA_C_BCGP)
    bcgp_pt(r->order, tt, &tp, &dtpdtt, &d2tpdtt2);
  else{
    tp = tt; dtpdtt = 1; d2tpdtt2 = 0;
  }

  if(p->info->number == XC_GGA_C_XPBE)
    gamm = cnst_beta*cnst_beta/(2.0*0.197363); /* for xPBE */
  else
    gamm = (1.0 - LOG(2.0))/(M_PI*M_PI);

  if(p->info->number == XC_GGA_C_PBELOC || p->info->number == XC_GGA_C_SG4) {
    /* PBEloc: \beta = \beta_0 + a t^2 f(r_s) */
    beta_den = EXP(- r->rs * r->rs);
    if(p->info->number == XC_GGA_C_PBELOC)
      beta = pbeloc_b0 + pbeloc_a * tt*tt * (1.0 - beta_den);
    else
      beta = sg4_b0 + sg4_sigma * tt * (1.0 - beta_den);

  }else if(p->info->number == XC_GGA_C_REGTPSS || p->func == 99)
    XC(beta_Hu_Langreth) (r->rs, r->order, &beta, &dbetadrs, &d2betadrs2);

  else
    beta = cnst_beta;

  if(p->info->number != XC_GGA_C_BCGP)
    pbe_eq8(r->order, beta, gamm, pw.zk, phi,
            &A, &dAdbeta, &dAdec, &dAdphi, &d2Adec2, &d2Adecphi, &d2Adphi2);
  else{ /* for BCGP */
    A = 1.0;
    dAdbeta = dAdec = dAdphi = d2Adec2 = d2Adecphi = d2Adphi2 = 0.0;
  }

  /* the sPBE functional contains one term less than the original PBE, so we set it to zero */
  B = (p->info->number == XC_GGA_C_SPBE) ? 0.0 : 1.0;
  pbe_eq7(r->order, p->func, beta, gamm, phi, tp, A, B,
	  &H, &dHdbeta, &dHdphi, &dHdt, &dHdA, &d2Hdphi2, &d2Hdphit, &d2HdphiA, &d2Hdt2, &d2HdtA, &d2HdA2);

  r->f = pw.zk + H;

  if(r->order < 1) return;

  /* full derivatives of functional */
  dfdbeta= dHdbeta + dHdA*dAdbeta;
  dfdphi = dHdphi + dHdA*dAdphi;
  dfdt   = dHdt*dtpdtt;
  dfdec  = 1.0 + dHdA*dAdec;

  dphidz = 0.0;
  if(auxp > p->info->min_zeta) dphidz += 1/auxp;
  if(auxm > p->info->min_zeta) dphidz -= 1/auxm;
  dphidz *= 1.0/3.0;

  dtdrs  = -r->xt/(2.0*tconv*phi*r->rs*pw.rs[0]);
  dtdxt  =  tt/r->xt;
  dtdphi = -tt/phi;

  if(p->info->number == XC_GGA_C_PBELOC){
    dfdt += dfdbeta * 2.0*pbeloc_a*tt*(1.0 - beta_den);
    dbetadrs = 2.0*pbeloc_a*tt*tt*r->rs*beta_den;

  }else if(p->info->number == XC_GGA_C_SG4){
    dfdt += dfdbeta * sg4_sigma*(1.0 - beta_den);
    dbetadrs = 2.0*sg4_sigma*tt*r->rs*beta_den;

  }else if(p->info->number != XC_GGA_C_REGTPSS && p->func != 99){
    dbetadrs = 0.0;
  }

  r->dfdrs   = dfdec*pw.dedrs + dfdt*dtdrs + dfdbeta*dbetadrs;
  r->dfdz    = dfdec*pw.dedz + (dfdphi + dfdt*dtdphi)*dphidz;
  r->dfdxt   = dfdt*dtdxt;
  r->dfdxs[0] = 0.0;
  r->dfdxs[1] = 0.0;

  if(r->order < 2) return;

  /* full derivatives of functional with respect to phi and zk */
  d2fdphi2  = d2Hdphi2 + 2.0*d2HdphiA*dAdphi + dHdA*d2Adphi2 + d2HdA2*dAdphi*dAdphi;
  d2fdphit  = (d2Hdphit + d2HdtA*dAdphi)*dtpdtt;
  d2fdphiec = d2HdphiA*dAdec + d2HdA2*dAdphi*dAdec + dHdA*d2Adecphi;
  d2fdt2    = d2Hdt2*dtpdtt*dtpdtt + dHdt*d2tpdtt2;
  d2fdtec   = d2HdtA*dAdec*dtpdtt;
  d2fdec2   = d2HdA2*dAdec*dAdec + dHdA*d2Adec2;

  d2phidz2 = 0.0;
  if(auxp > p->info->min_zeta) d2phidz2 += 1.0/((1.0 + r->zeta)*auxp);
  if(auxm > p->info->min_zeta) d2phidz2 += 1.0/((1.0 - r->zeta)*auxm);
  d2phidz2 *= -1.0/9.0;

  d2tdrs2   =  3.0*r->xt/(4.0*tconv*phi*pw.rs[2]*pw.rs[0]);
  d2tdrsxt  =  dtdrs/r->xt;
  d2tdphi2  = -2.0*dtdphi/phi;
  d2tdrsphi = -dtdrs/phi;
  d2tdxtphi =  dtdphi/r->xt;

  r->d2fdrs2     = dfdec*pw.d2edrs2 + d2fdec2*pw.dedrs*pw.dedrs + 2.0*d2fdtec*pw.dedrs*dtdrs + d2fdt2*dtdrs*dtdrs + dfdt*d2tdrs2;
  r->d2fdrsz     = dfdec*pw.d2edrsz + pw.dedrs*(d2fdec2*pw.dedz + dphidz*(d2fdtec*dtdphi + d2fdphiec))
    + dfdt*dphidz*d2tdrsphi + dtdrs*(d2fdtec*pw.dedz + dphidz*(d2fdt2*dtdphi + d2fdphit));
  r->d2fdrsxt    = dtdxt*(d2fdtec*pw.dedrs + d2fdt2*dtdrs) + dfdt*d2tdrsxt;
  r->d2fdrsxs[0] = 0.0;
  r->d2fdrsxs[1] = 0.0;
  r->d2fdz2      = dfdec*pw.d2edz2 + d2fdec2*pw.dedz*pw.dedz + dfdt*(dtdphi*d2phidz2 + d2tdphi2*dphidz*dphidz)
    + dfdphi*d2phidz2 + 2.0*dphidz*pw.dedz*(d2fdtec*dtdphi + d2fdphiec)
    + dphidz*dphidz*(d2fdt2*dtdphi*dtdphi + 2.0*d2fdphit*dtdphi + d2fdphi2);
  r->d2fdzxt     = dfdt*d2tdxtphi*dphidz + dtdxt*(d2fdtec*pw.dedz + dphidz*(d2fdt2*dtdphi + d2fdphit));
  r->d2fdzxs[0]  = 0.0;
  r->d2fdzxs[1]  = 0.0;
  r->d2fdxt2     = d2fdt2*dtdxt*dtdxt;
  r->d2fdxtxs[0] = 0.0;
  r->d2fdxtxs[1] = 0.0;
  r->d2fdxs2[0]  = 0.0;
  r->d2fdxs2[1]  = 0.0;
  r->d2fdxs2[2]  = 0.0;
}

#define func XC(gga_c_pbe_func)
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_pbe) = {
  XC_GGA_C_PBE,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1996_3865, &xc_ref_Perdew1996_3865_err, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};


const XC(func_info_type) XC(func_info_gga_c_pbe_sol) = {
  XC_GGA_C_PBE_SOL,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew2008_136406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};


const XC(func_info_type) XC(func_info_gga_c_xpbe) = {
  XC_GGA_C_XPBE,
  XC_CORRELATION,
  "Extended PBE by Xu & Goddard III",
  XC_FAMILY_GGA,
  {&xc_ref_Xu2004_4068, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_pbe_jrgx) = {
  XC_GGA_C_PBE_JRGX,
  XC_CORRELATION,
  "Reparametrized PBE by Pedroza, Silva & Capelle",
  XC_FAMILY_GGA,
  {&xc_ref_Pedroza2009_201106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_rge2) = {
  XC_GGA_C_RGE2,
  XC_CORRELATION,
  "Regularized PBE",
  XC_FAMILY_GGA,
  {&xc_ref_Ruzsinszky2009_763, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_apbe) = {
  XC_GGA_C_APBE,
  XC_CORRELATION,
  "mu fixed from the semiclassical neutral atom",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_186406, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_spbe) = {
  XC_GGA_C_SPBE,
  XC_CORRELATION,
  "PBE correlation to be used with the SSB exchange",
  XC_FAMILY_GGA,
  {&xc_ref_Swart2009_094103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_regtpss) = {
  XC_GGA_C_REGTPSS,
  XC_CORRELATION,
  "regularized TPSS correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew2009_026403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_zpbesol) = {
  XC_GGA_C_ZPBESOL,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEsol",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};


const XC(func_info_type) XC(func_info_gga_c_pbeint) = {
  XC_GGA_C_PBEINT,
  XC_CORRELATION,
  "PBE for hybrid interfaces",
  XC_FAMILY_GGA,
  {&xc_ref_Fabiano2010_113104, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_zpbeint) = {
  XC_GGA_C_ZPBEINT,
  XC_CORRELATION,
  "spin-dependent gradient correction to PBEint",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2011_233103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_pbeloc) = {
  XC_GGA_C_PBELOC,
  XC_CORRELATION,
  "Semilocal dynamical correlation",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2012_035130, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_bcgp) = {
  XC_GGA_C_BCGP,
  XC_CORRELATION,
  "Burke, Cancio, Gould, and Pittalis",
  XC_FAMILY_GGA,
  {&xc_ref_Burke2014_4834, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-12, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_pbefe) = {
  XC_GGA_C_PBEFE,
  XC_CORRELATION,
  "PBE for formation energies",
  XC_FAMILY_GGA,
  {&xc_ref_Perez2015_3844, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_c_pbe_mol) = {
  XC_GGA_C_PBE_MOL,
  XC_EXCHANGE,
  "Reparametrized PBE by del Campo, Gazquez, Trickey & Vela",
  XC_FAMILY_GGA,
  {&xc_ref_delCampo2012_104108, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, 
  NULL, work_gga_c, NULL
};


const XC(func_info_type) XC(func_info_gga_c_sg4) = {
  XC_GGA_C_SG4,
  XC_EXCHANGE,
  "Semiclassical GGA at fourth order",
  XC_FAMILY_GGA,
  {&xc_ref_Constantin2016_045126, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_c_pbe_init, NULL, NULL,
  work_gga_c,
  NULL
};
