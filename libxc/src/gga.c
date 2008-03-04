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
#include "funcs_gga.c"


/* initialization */
int XC(gga_init)(XC(gga_type) *p, int functional, int nspin)
{
  int i;

  assert(p != NULL);

  /* let us first find out if we know the functional */
  for(i=0; XC(gga_known_funct)[i]!=NULL; i++){
    if(XC(gga_known_funct)[i]->number == functional) break;
  }
  if(XC(gga_known_funct)[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->params = NULL;
  p->info = XC(gga_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  /* see if we need to initialize the functional */
  if(p->info->init != NULL)
    p->info->init(p);
  return 0;
}


/* Termination */
void XC(gga_end)(XC(gga_type) *p)
{
  assert(p != NULL);

  if(p->info->end != NULL)
    p->info->end(p);
}


/* Some useful formulas:

   sigma_st  = grad rho_s . grad rho_t
   vrho_s    = d e / d rho_s
   vsigma_st = d n*e / d sigma_st

if nspin == 2
   rho   = (rho_u, rho_d)
   sigma = (sigma_uu, sigma_du, sigma_dd)
*/
void XC(gga)(XC(gga_type) *p, FLOAT *rho, FLOAT *sigma,
	    FLOAT *e, FLOAT *vrho, FLOAT *vsigma)
{
  FLOAT dens;
  int i;

  assert(p!=NULL);
  
  *e = 0.0;
  for(i=0; i<p->nspin; i++) vrho [i] = 0.0;

  vsigma[0] = 0.0;
  if(p->nspin == XC_POLARIZED){
    vsigma[1] = 0.0; vsigma[2] = 0.0;
  }

  dens = rho[0];
  if(p->nspin == XC_POLARIZED) dens += rho[1];
  if(dens <= MIN_DENS) return;

  assert(p->info!=NULL && p->info->gga!=NULL);
  p->info->gga(p, rho, sigma, e, vrho, vsigma);
}

