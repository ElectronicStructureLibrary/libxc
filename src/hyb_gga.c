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
#include "funcs_hyb_gga.c"

/* initialization */
/*****************************************************/
int XC(hyb_gga_init)(XC(hyb_gga_type) *p, int functional, int nspin)
{
  int i;

  assert(p!=NULL);

  /* let us first find out if we know the functional */
  for(i=0; XC(hyb_gga_known_funct)[i]!=NULL; i++){
    if(XC(hyb_gga_known_funct)[i]->number == functional) break;
  }
  if(XC(hyb_gga_known_funct)[i] == NULL) return -1; /* functional not found */

  /* initialize structure */
  p->params = NULL;
  p->mix    = NULL;
  p->info = XC(hyb_gga_known_funct)[i];

  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;

  p->exx_coef = 1.0;

  /* we always need to initialize the functional */
  assert(p->info->init != NULL);
  p->info->init(p);
  return 0;
}


/* Termination */
/*****************************************************/
void XC(hyb_gga_end)(XC(hyb_gga_type) *p)
{
  assert(p!=NULL);

  if(p->info->end != NULL)
    p->info->end(p);

  if(p->mix != NULL)
    XC(mix_func_free)(p->mix);
  free(p->mix); p->mix = NULL;
}


/*****************************************************/
void XC(hyb_gga)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  assert(p!=NULL && p->info!=NULL);
  
  if(!XC(gga_input_init)(p->info, p->nspin, rho, zk, vrho, vsigma,
			 v2rho2, v2rhosigma, v2sigma2)) return; 

  /* hybrid may want to add some term */
  if(p->info->gga!=NULL)
    p->info->gga(p, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

  /* and now we have the general mixture of GGAs */
  if(p->mix != NULL)
    XC(mix_func)(p->mix, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}


/* especializations */
inline void 
XC(hyb_gga_exc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		FLOAT *zk)
{
  XC(hyb_gga)(p, rho, sigma, zk, NULL, NULL, NULL, NULL, NULL);
}

inline void 
XC(hyb_gga_exc_vxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(hyb_gga)(p, rho, sigma, zk, vrho, vsigma, NULL, NULL, NULL);
}

inline void 
XC(hyb_gga_vxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *vrho, FLOAT *vsigma)
{
  XC(hyb_gga)(p, rho, sigma, NULL, vrho, vsigma, NULL, NULL, NULL);
}

inline void 
XC(hyb_gga_fxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(hyb_gga)(p, rho, sigma, NULL, NULL, NULL, v2rho2, v2rhosigma, v2sigma2);
}


/*****************************************************/
FLOAT XC(hyb_gga_exx_coef)(XC(hyb_gga_type) *p)
{
  assert(p!=NULL);

  return p->exx_coef;
}
