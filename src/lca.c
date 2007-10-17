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

#define XC_LCA_OMC       301   /* Orestes, Marcasso & Capelle  */
#define XC_LCA_LCH       302   /* Lee, Colwell & Handy         */

/* initialization */
void xc_lca_init(xc_lca_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_LCA_LCH ||
	 functional == XC_LCA_OMC);

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals */
  switch(functional){
  case XC_LCA_LCH:
    lca_lch_init(p);
    break;

  case XC_LCA_OMC:
    lca_omc_init(p);
    break;
  }

}

/* WARNING - should use new definition of input/output !! */
#define   _(is, x)   [3*is + x]

void xc_lca(xc_lca_type *p, double *rho, double *v, double *e, double *dedd, double *dedv)
{
  int i;
  int j;
  double rs, drsdd, vs, vs2, kf, dkfdrs, f, s, dsdrs;
  double dens;

  assert(p!=NULL);

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  
  if(dens <= MIN_DENS){
    *e = 0.0;
    for(i=0; i<  p->nspin; i++) dedd [i] = 0.0;
    for(i=0; i<3*p->nspin; i++) dedv [i] = 0.0;
    return;
  }

  for(i=0; i<  p->nspin; i++){
    
    rs = RS(rho[i]);
    drsdd = -rs/(3.0*rho[i]);
    vs = sqrt(v _(i, 0)*v _(i, 0) + 
	      v _(i, 1)*v _(i, 1) +
	      v _(i, 2)*v _(i, 2));
    vs2 = vs*vs;
    dkfdrs = pow(9.0*M_PI/4.0, 1.0/3.0);
    kf = dkfdrs*rs;
    f = 24.0*M_PI*M_PI;

    switch(p->info->number){
    case XC_LCA_LCH:
      lca_s_lch(rs, &s, &dsdrs);
      break;

    case XC_LCA_OMC:
      lca_s_omc(rs, &s, &dsdrs);
      break;
    }

    *e = kf/f*(s - 1.0)*vs2;
    dedd [i] = drsdd/f*( dkfdrs*(s - 1.0) + kf*dsdrs )*vs2;
    for(j=0; j<3; j++) dedv _(i, j) = 2.0*kf/f*(s - 1.0)*v _(i, j);
  }
  
}

void xc_lca_sp(xc_lca_type *p, float *rho, float *v, 
	       float *e, float *dedd, float *dedv)
{
  double drho[2], dv[6];
  double de[1], ddedd[2], ddedv[6];
  int ii;

  for(ii=0; ii < p->nspin; ii++) drho[ii] = rho[ii];
  for(ii=0; ii < 3*p->nspin; ii++) dv[ii] = v[ii];
  
  xc_lca(p, drho, dv, de, ddedd, ddedv);

  e[0] = de[0];
  for(ii=0; ii < p->nspin; ii++) dedd[ii] = ddedd[ii];
  for(ii=0; ii < 3*p->nspin; ii++) dedv[ii] = ddedv[ii];
}
