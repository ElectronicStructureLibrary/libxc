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

/************************************************************************
 Implements Perdew, Burke & Ernzerhof Generalized Gradient Approximation
 correlation functional.

 I based this implementation on a routine from L.C. Balbas and J.M. Soler
************************************************************************/

#define XC_GGA_C_PBE          130 /* Perdew, Burke & Ernzerhof correlation     */
#define XC_GGA_C_PBE_SOL      133 /* Perdew, Burke & Ernzerhof correlation SOL */

static const FLOAT beta[2]  = {
  0.06672455060314922,  /* original PBE */
  0.046                 /* PBE sol      */
};
static const FLOAT gamm  = 0.03109069086965489503494086371273; /* (1.0 - log(2.0))/(M_PI*M_PI) */

static void gga_c_pbe_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_init)(p->lda_aux, XC_LDA_C_PW_MOD, p->nspin);
}

static void gga_c_pbe_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  free(p->lda_aux);
}


static inline void pbe_eq8(int func, FLOAT ecunif, FLOAT phi, 
		    FLOAT *A, FLOAT *dec, FLOAT *dphi)
{
  FLOAT phi3, f1, f2, f3, dx;
  static const FLOAT beta[2]  = {
    0.06672455060314922,  /* original PBE */
    0.046                 /* PBE sol      */
  };
  static const FLOAT gamm  = 0.03109069086965489503494086371273; /* (1.0 - log(2.0))/(M_PI*M_PI) */

  phi3 = POW(phi, 3);
  f1   = ecunif/(gamm*phi3);
  f2   = exp(-f1);
  f3   = f2 - 1.0;

  *A   = beta[func]/(gamm*f3);

  dx    = beta[func]*f2/(gamm*f3*f3);
  *dec  =  dx/(gamm*phi3);
  *dphi = -dx*3.0*ecunif/(gamm*phi*phi3);
}


static inline void pbe_eq7(int func, FLOAT phi, FLOAT t, FLOAT A, 
		    FLOAT *H, FLOAT *dphi, FLOAT *dt, FLOAT *dA)
{
  FLOAT t2, phi3, f1, f2, f3;
  static const FLOAT beta[2]  = {
    0.06672455060314922,  /* original PBE */
    0.046                 /* PBE sol      */
  };
  static const FLOAT gamm  = 0.03109069086965489503494086371273; /* (1.0 - log(2.0))/(M_PI*M_PI) */

  t2   = t*t;
  phi3 = POW(phi, 3);

  f1 = t2 + A*t2*t2;
  f3 = 1.0 + A*f1;
  f2 = beta[func]*f1/(gamm*f3);

  *H = gamm*phi3*log(1.0 + f2);

  {
    FLOAT df1dt, df2dt, df1dA, df2dA;

    *dphi = 3.0*gamm*phi*phi*log(1.0 + f2);
    
    df1dt = t*(2.0 + 4.0*A*t2);
    df2dt = (beta[func]/gamm) / (f3*f3) * df1dt;
    *dt   = gamm*phi3*df2dt/(1.0 + f2);
    
    df1dA = t2*t2;
    df2dA = (beta[func]/gamm) / (f3*f3) * (df1dA - f1*f1);
    *dA   = gamm*phi3*df2dA/(1.0 + f2);
  }

}

static void gga_c_pbe(void *p_, FLOAT *rho, FLOAT *sigma,
		      FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;
  XC(perdew_t) pt;

  int func;
  FLOAT A, dAdec, dAdphi;
  FLOAT H, dHdphi, dHdt, dHdA;

  func = (p->info->number == XC_GGA_C_PBE_SOL) ? 1 : 0;

  XC(perdew_params)(p, rho, sigma, &pt);

  pbe_eq8(func, pt.ecunif, pt.phi, &A, &dAdec, &dAdphi);
  pbe_eq7(func, pt.phi, pt.t, A, &H, &dHdphi, &dHdt, &dHdA);

  *e = pt.ecunif + H;

  pt.dphi    = dHdphi + dHdA*dAdphi;
  pt.dt      = dHdt;
  pt.decunif = 1.0 + dHdA*dAdec;

  XC(perdew_potentials)(&pt, rho, *e, vrho, vsigma);
}


const XC(func_info_type) XC(func_info_gga_c_pbe) = {
  XC_GGA_C_PBE,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof",
  XC_FAMILY_GGA,
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_c_pbe_init,
  gga_c_pbe_end,
  NULL,            /* this is not an LDA                   */
  gga_c_pbe,
};

const XC(func_info_type) XC(func_info_gga_c_pbe_sol) = {
  XC_GGA_C_PBE_SOL,
  XC_CORRELATION,
  "Perdew, Burke & Ernzerhof SOL",
  XC_FAMILY_GGA,
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)\n"
  "JP Perdew, K Burke, and M Ernzerhof, Phys. Rev. Lett. 78, 1396(E) (1997)"
  "JP Perdew, et al, arXiv:0707.2088v1",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_c_pbe_init,
  gga_c_pbe_end,
  NULL,            /* this is not an LDA                   */
  gga_c_pbe,
};
