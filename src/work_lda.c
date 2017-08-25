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

/**
 * @file work_lda.c
 * @brief This file is to be included in LDA functionals. As often these
 *        functionals are written as a function of rs and zeta, this
 *        routine performs the necessary conversions between this and a functional
 *        of rho.
 */


#ifndef XC_DIMENSIONS
#define XC_DIMENSIONS 3
#endif


/**
 * @param[in,out] func_type: pointer to pspdata structure to be initialized
 */
static void 
work_lda(const XC(func_type) *p, int np, const double *rho, 
	 double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  XC(lda_work_t) r;
  int is, ip;
  double dens, drs, d2rs, d3rs;

  /* Wigner radius */
# if   XC_DIMENSIONS == 1
  const double cnst_rs = 0.5;
# elif XC_DIMENSIONS == 2
  const double cnst_rs = 1.0/M_SQRTPI;
# else /* three dimensions */
  const double cnst_rs = RS_FACTOR;
# endif

  /* Initialize memory */
  memset(&r, 0, sizeof(r));

  r.order = -1;
  if(zk     != NULL) r.order = 0;
  if(vrho   != NULL) r.order = 1;
  if(v2rho2 != NULL) r.order = 2;
  if(v3rho3 != NULL) r.order = 3;
  if(r.order < 0) return;

  for(ip = 0; ip < np; ip++){
    XC(rho2dzeta)(p->nspin, rho, &dens, &r.z);

    if(dens < p->dens_threshold) goto end_ip_loop;

    r.rs = cnst_rs*POW(dens, -1.0/XC_DIMENSIONS);

    func(p, &r);

    if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
      *zk = r.f;

    if(r.order < 1) goto end_ip_loop;

    drs = -r.rs/(XC_DIMENSIONS*dens);
    
    if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
      vrho[0] = r.f + dens*r.dfdrs*drs;

      if(p->nspin == XC_POLARIZED){
	vrho[1] = vrho[0] - (r.z + 1.0)*r.dfdz;
	vrho[0] = vrho[0] - (r.z - 1.0)*r.dfdz;
      }
    }
  
    if(r.order < 2) goto end_ip_loop;
    
    d2rs = -drs*(1.0 + XC_DIMENSIONS)/(XC_DIMENSIONS*dens);
    
    if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
      v2rho2[0] = r.dfdrs*(2.0*drs + dens*d2rs) + dens*r.d2fdrs2*drs*drs;
      
      if(p->nspin == XC_POLARIZED){
	double sign[3][2] = {{-1.0, -1.0}, {-1.0, +1.0}, {+1.0, +1.0}};
	
	for(is=2; is>=0; is--){
	  v2rho2[is] = v2rho2[0] - r.d2fdrsz*(2.0*r.z + sign[is][0] + sign[is][1])*drs
	    + (r.z + sign[is][0])*(r.z + sign[is][1])*r.d2fdz2/dens;
	}
      }
    }
    
    if(r.order < 3) goto end_ip_loop;

    d3rs = -d2rs*(1.0 + 2.0*XC_DIMENSIONS)/(XC_DIMENSIONS*dens);
    
    if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC)){
      v3rho3[0] = r.dfdrs*(3.0*d2rs + dens*d3rs) + 
	3.0*r.d2fdrs2*drs*(drs + dens*d2rs) + r.d3fdrs3*dens*drs*drs*drs;
      
      if(p->nspin == XC_POLARIZED){
	double sign[4][3] = {{-1.0, -1.0, -1.0}, {-1.0, -1.0, +1.0}, {-1.0, +1.0, +1.0}, {+1.0, +1.0, +1.0}};
	
	for(is=3; is>=0; is--){
	  double ff;
	  
	  v3rho3[is]  = v3rho3[0] - (2.0*r.z  + sign[is][0] + sign[is][1])*(d2rs*r.d2fdrsz + drs*drs*r.d3fdrs2z);
	  v3rho3[is] += (r.z + sign[is][0])*(r.z + sign[is][1])*(-r.d2fdz2/dens + r.d3fdrsz2*drs)/dens;
	  
	  ff  = r.d2fdrsz*(2.0*drs + dens*d2rs) + dens*r.d3fdrs2z*drs*drs;
	  ff += -2.0*r.d2fdrsz*drs - r.d3fdrsz2*(2.0*r.z + sign[is][0] + sign[is][1])*drs;
	  ff += (r.z + sign[is][0])*(r.z + sign[is][1])*r.d3fdz3/dens;
	  ff += (2.0*r.z  + sign[is][0] + sign[is][1])*r.d2fdz2/dens;
	  
	  v3rho3[is] += -ff*(r.z + sign[is][2])/dens;
	}
      }
    }

  end_ip_loop:
    rho += p->n_rho;

    if(zk != NULL)
      zk += p->n_zk;
    
    if(vrho != NULL)
      vrho += p->n_vrho;

    if(v2rho2 != NULL)
      v2rho2 += p->n_v2rho2;

    if(v3rho3 != NULL)
      v3rho3 += p->n_v3rho3;

  } /* for(ip) */
}
