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
Correlation functional by Pittalis, Rasanen & Marques for the 2D electron gas
************************************************************************/

#define XC_LDA_C_2D_PRM  16   /* Pittalis, Rasanen & Marques correlation in 2D */

typedef struct{
  FLOAT N;
  FLOAT c;
} lda_c_prm_params;

/* parameters necessary to the calculation */
static FLOAT prm_q = 3.9274; /* 2.258 */

/* Initialization */
static void
lda_c_2d_prm_init(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;
  lda_c_prm_params *params;

  assert(p->params == NULL);

  p->params = malloc(sizeof(lda_c_prm_params));
  params = (lda_c_prm_params *) (p->params);

  params->N = 0.0;
}


static void 
lda_c_2d_prm_end(void *p_)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;

  assert(p->params != NULL);
  free(p->params);
  p->params = NULL;
}


void 
XC(lda_c_2d_prm_set_params)(XC(lda_type) *p, FLOAT N)
{
  lda_c_prm_params *params;

  assert(p->params != NULL);
  params = (lda_c_prm_params *) (p->params);

  if(N <= 1){
    fprintf(stderr, "PRM functional can not be used for N_electrons <= 1\n");
    exit(1);
  }

  params->N = N;
  params->c = M_PI/(2.0*(N - 1.0)*prm_q*prm_q); /* Eq. (13) */
}


static void
lda_c_2d_prm(const void *p_, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2)
{
  XC(lda_type) *p = (XC(lda_type) *)p_;
  lda_c_prm_params *params;

  FLOAT dens, beta, phi, c;
  FLOAT sqpi, t1, t2, t3, dt1dbeta, dt1dphi, dt3dphi, dbetadn, dphidn;

  assert(p != NULL && p->params != NULL);
  params = p->params;

  assert(params->N > 1.0);
  assert(p->nspin == XC_UNPOLARIZED);
  
  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];

  beta = prm_q*sqrt(dens); /* Eq. (4) */

  sqpi = sqrt(M_PI);
  c    = params->c;

  phi = beta/(beta + sqpi/2.0);

  t3  = phi - 1.0; /* original version has (phi-1)^2 */
  t2  = M_PI/(2.0*prm_q*prm_q);

  t1  = sqpi*beta*t3/(2.0*sqrt(2.0 + c));
  t1 += phi*(phi - 1.0)/(2.0 + c);
  t1 += sqpi*phi*phi/(4.0*beta*POW(2.0 + c, 1.5));
  t1 += sqpi*beta*(phi - 1.0)/sqrt(1.0 + c);
  t1 += phi/(1.0 + c);
  t1 *= t2;

  if(zk != NULL)
    *zk = t1;

  if(vrho == NULL) return;

  dt1dbeta  = sqpi*t3/(2.0*sqrt(2.0 + c));
  dt1dbeta -= sqpi*phi*phi/(4.0*beta*beta*POW(2.0 + c, 1.5));
  dt1dbeta += sqpi*(phi - 1.0)/sqrt(1.0 + c);
  dt1dbeta *= t2;

  dt3dphi   = 1.0;
  dt1dphi   = sqpi*beta/(2.0*sqrt(2.0 + c))*dt3dphi;
  dt1dphi  += (2.0*phi - 1.0)/(2.0 + c);
  dt1dphi  += sqpi*2.0*phi/(4.0*beta*POW(2.0 + c, 1.5));
  dt1dphi  += sqpi*beta/sqrt(1.0 + c);
  dt1dphi  += 1.0/(1.0 + c);
  dt1dphi  *= t2;

  dbetadn   = prm_q/(2.0*sqrt(dens));
  dphidn    = sqpi/(2.0*(beta + sqpi/2.0)*(beta + sqpi/2.0));
  dphidn   *= dbetadn;

  *vrho     = t1 + dens*(dt1dbeta*dbetadn + dt1dphi*dphidn);
}


const XC(func_info_type) XC(func_info_lda_c_2d_prm) = {
  XC_LDA_C_2D_PRM,
  XC_CORRELATION,
  "PRM (for 2D systems)",
  XC_FAMILY_LDA,
  "S Pittalis, E Rasanen, and MAL Marques, Phys. Rev. B 78, 195322 (2008)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  lda_c_2d_prm_init,
  lda_c_2d_prm_end,
  lda_c_2d_prm
};
