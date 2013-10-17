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

/************************************************************************
  This file is to be included in GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3), this
  routine performs the necessary conversions between a functional of s
  and of rho.
************************************************************************/

#ifndef HEADER
#  define HEADER 1
#endif

#ifndef XC_DIMENSIONS
#  define XC_DIMENSIONS 3
#endif

static void
#ifdef XC_KINETIC_FUNCTIONAL
work_gga_k
#else
work_gga_x
#endif
(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
 FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
{
  FLOAT sfact, sfact2, x_factor_c, alpha, beta, dens;
  int is, is2, ip, order;

  /* constants for the evaluation of the different terms */
  FLOAT c_zk[1];
  FLOAT c_vrho[3], c_vsigma[2];
  FLOAT c_v2rho2[3], c_v2rhosigma[4], c_v2sigma2[2];
  FLOAT c_v3rho3[4], c_v3rho2sigma[3], c_v3rhosigma2[3], c_v3sigma3[3];

  /* variables used inside the is loop */
  FLOAT gdm, ds, rhoLDA;
  FLOAT x, f, dfdx, d2fdx2, d3fdx3, lvsigma, lv2sigma2, lvsigmax, lvrho;

  /* alpha is the power of rho in the corresponding LDA
     beta  is the power of rho in the expression for x */

  beta = 1.0 + 1.0/XC_DIMENSIONS; /* exponent of the density in expression for x */

#ifndef XC_KINETIC_FUNCTIONAL
  alpha = beta;

#  if XC_DIMENSIONS == 2
  x_factor_c = -X_FACTOR_2D_C;
#  else /* three dimensions */
  x_factor_c = -X_FACTOR_C;
#  endif

#else

#  if XC_DIMENSIONS == 2
#  else /* three dimensions */
  alpha = 5.0/3.0;
  x_factor_c = K_FACTOR_C;
#  endif

#endif

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  /* Initialize several constants */
  order = -1;
  if(zk     != NULL){
    order = 0;
    c_zk[0] = sfact*x_factor_c;
  }
  if(vrho   != NULL){
    order = 1;
    c_vrho[0]   =  x_factor_c*alpha;
    c_vrho[1]   = -x_factor_c*beta;
    c_vrho[2]   =  x_factor_c;
    c_vsigma[0] =  sfact*x_factor_c;
    c_vsigma[1] =  sfact*x_factor_c;
  }
  if(v2rho2 != NULL){
    order = 2;
    c_v2rho2[0] = (x_factor_c/sfact) * (alpha - 1.0)*alpha;
    c_v2rho2[1] = (x_factor_c/sfact) * beta*(beta - 2.0*alpha + 1.0);
    c_v2rho2[2] = (x_factor_c/sfact) * beta*beta;
    c_v2rhosigma[0] =  x_factor_c * (alpha - beta)/2.0;
    c_v2rhosigma[1] = -x_factor_c * beta/2.0;
    c_v2rhosigma[2] =  x_factor_c * alpha;
    c_v2rhosigma[3] = -x_factor_c * beta;
    c_v2sigma2[0] = x_factor_c*sfact / 4.0; 
    c_v2sigma2[1] = x_factor_c*sfact;
  }
  if(v3rho3 != NULL){
    order = 3;
    c_v3rho3[0] =  (x_factor_c/sfact2) * (alpha - 2.0)*(alpha - 1.0)*alpha;
    c_v3rho3[1] = -(x_factor_c/sfact2) * (3.0*alpha*alpha - 3.0*alpha*(2.0 + beta) + (1.0 + beta)*(2.0 + beta))*beta;
    c_v3rho3[2] = -(x_factor_c/sfact2) * 3.0*(1.0 - alpha + beta)*beta*beta;
    c_v3rho3[3] = -(x_factor_c/sfact2) * beta*beta*beta;
    c_v3rho2sigma[0] = (x_factor_c/sfact) * (alpha - beta - 1.0)*(alpha - beta)/2.0;
    c_v3rho2sigma[1] = (x_factor_c/sfact) * (1.0 - 2.0*alpha + 3.0*beta)*beta/2.0;
    c_v3rho2sigma[2] = (x_factor_c/sfact) * beta*beta/2.0;
    c_v3rhosigma2[0] = -x_factor_c * (alpha - beta)/4.0;
    c_v3rhosigma2[1] =  x_factor_c * (alpha - beta)/4.0;
    c_v3rhosigma2[2] = -x_factor_c * beta/4.0;
    c_v3sigma3[0] =  x_factor_c*sfact * 3.0/8.0;
    c_v3sigma3[1] = -x_factor_c*sfact * 3.0/8.0;
    c_v3sigma3[2] =  x_factor_c*sfact /8.0;
  }
  if(order < 0) return;


  /* the loop over the points starts */
  for(ip = 0; ip < np; ip++){
    dens = (p->nspin == XC_UNPOLARIZED) ? rho[0] : rho[0] + rho[1];
    if(dens < p->info->min_dens) goto end_ip_loop;

    for(is=0; is<p->nspin; is++){
      is2 = 2*is;

      if(rho[is] < p->info->min_dens) continue;

      gdm    = max(SQRT(sigma[is2])/sfact, p->info->min_grad);
      ds     = rho[is]/sfact;
      rhoLDA = POW(ds, alpha);
      x      = gdm/POW(ds, beta);

      dfdx = d2fdx2 = d3fdx3 = 0.0;
      lvsigma = lv2sigma2 = lvsigmax = lvrho = 0.0;

#if   HEADER == 1

      func(p, order, x, &f, &dfdx, &d2fdx2, &d3fdx3);

#elif HEADER == 2

      /* this second header is useful for functionals that depend
	 explicitly both on x and on sigma */
      func(p, order, x, gdm*gdm, &f, &dfdx, &lvsigma, &d2fdx2, &lv2sigma2, &lvsigmax);
      
      lvsigma   /= sfact2;
      lvsigmax  /= sfact2;
      lv2sigma2 /= sfact2*sfact2;

#elif HEADER == 3

      /* this second header is useful for functionals that depend
	 explicitly both on x and on rho*/
      func(p, order, x, ds, &f, &dfdx, &lvrho);

#endif

      if(order > 0) dfdx   *= x;
      if(order > 1) d2fdx2 *= x*x;
      if(order > 2) d3fdx3 *= x*x*x;

      if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
	*zk += rhoLDA*
	  c_zk[0]*f;
      
      if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
	vrho[is] += (rhoLDA/ds)*
	  (c_vrho[0]*f + c_vrho[1]*dfdx) + rhoLDA*c_vrho[2]*lvrho;
	
	if(gdm > p->info->min_grad)
	  vsigma[is2] = rhoLDA*
	    (c_vsigma[0]*dfdx/(2.0*sigma[is2]) + c_vsigma[1]*lvsigma);
      }
      
      if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
	v2rho2[is2] = rhoLDA/(ds*ds) * (c_v2rho2[0]*f + c_v2rho2[1]*dfdx + c_v2rho2[2]*d2fdx2);
	
	if(gdm > p->info->min_grad){
	  v2rhosigma[is*5] = (rhoLDA/ds) *
	    ((c_v2rhosigma[0]*dfdx + c_v2rhosigma[1]*d2fdx2)/sigma[is2] + c_v2rhosigma[2]*lvsigma + c_v2rhosigma[3]*x*lvsigmax);
	  v2sigma2  [is*5] = rhoLDA*
	    (c_v2sigma2[0]*(d2fdx2 - dfdx)/(sigma[is2]*sigma[is2]) + c_v2sigma2[1]*(lv2sigma2 + lvsigmax*x/sigma[is2]));
	}
      }

      if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC)){
	v3rho3[is*3] = rhoLDA/(ds*ds*ds) *
	  (c_v3rho3[0]*f + c_v3rho3[1]*dfdx + c_v3rho3[2]*d2fdx2 + c_v3rho3[3]*d3fdx3);

	if(gdm > p->info->min_grad){
	  v3rho2sigma[is*8] = rhoLDA/(ds*ds) *
	    (c_v3rho2sigma[0]*dfdx + c_v3rho2sigma[1]*d2fdx2 + c_v3rho2sigma[2]*d3fdx3)/sigma[is2];

	  v3rhosigma2[is*11] = (rhoLDA/ds) *
	    (c_v3rhosigma2[0]*dfdx + c_v3rhosigma2[1]*d2fdx2 + c_v3rhosigma2[2]*d3fdx3)/(sigma[is2]*sigma[is2]);

	  v3sigma3[is*9] = rhoLDA*
	    (c_v3sigma3[0]*dfdx + c_v3sigma3[1]*d2fdx2 + c_v3sigma3[2]*d3fdx3)/(sigma[is2]*sigma[is2]*sigma[is2]);
	}
      }
    }

    if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
      *zk /= dens; /* we want energy per particle */
    
  end_ip_loop:
    /* increment pointers */
    rho   += p->n_rho;
    sigma += p->n_sigma;
    
    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL){
      vrho   += p->n_vrho;
      vsigma += p->n_vsigma;
    }

    if(v2rho2 != NULL){
      v2rho2     += p->n_v2rho2;
      v2rhosigma += p->n_v2rhosigma;
      v2sigma2   += p->n_v2sigma2;
    }

    if(v3rho3 != NULL){
      v3rho3      += p->n_v3rho3;
      v3rho2sigma += p->n_v3rho2sigma;
      v3rhosigma2 += p->n_v3rhosigma2;
      v3sigma3    += p->n_v3sigma3;
    }
  }
}
