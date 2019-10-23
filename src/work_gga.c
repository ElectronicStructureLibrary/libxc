/*
 Copyright (C) 2006-2018 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_gga.c
 * @brief This file is to be included in GGA functionals.
 */

/* hack to avoid compiler warnings */
#define NOARG

#ifdef XC_NO_EXC
#define OUT_PARAMS GGA_OUT_PARAMS_NO_EXC(NOARG)
#else
#define OUT_PARAMS zk, GGA_OUT_PARAMS_NO_EXC(NOARG)
#endif

/**
 * @param[in,out] func_type: pointer to functional structure
 */
static void 
work_gga(const XC(func_type) *p, int np,
         const double *rho, const double *sigma,
         double *zk, double *vrho, double *vsigma,
         double *v2rho2, double *v2rhosigma, double *v2sigma2,
         double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  int ip, order;
  double dens, zeta;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;

  if(order < 0) return;

  for(ip = 0; ip < np; ip++){
    xc_rho2dzeta(p->nspin, rho, &dens, &zeta);

    if(dens > p->dens_threshold){
      if(p->nspin == XC_UNPOLARIZED){             /* unpolarized case */
        func_unpol(p, order, rho, sigma, OUT_PARAMS);
      
      }else if(zeta >  1.0 - 1e-10){              /* ferromagnetic case - spin 0 */
        func_ferr(p, order, rho, sigma, OUT_PARAMS);
        
      }else if(zeta < -1.0 + 1e-10){              /* ferromagnetic case - spin 1 */
        internal_counters_gga_next
          (&(p->dim), -1, &rho, &sigma, &zk, &vrho, &vsigma,
           &v2rho2, &v2rhosigma, &v2sigma2, &v3rho3, &v3rho2sigma, &v3rhosigma2, &v3sigma3);
        func_ferr(p, order, rho, sigma, OUT_PARAMS);
        internal_counters_gga_prev
          (&(p->dim), -1, &rho, &sigma, &zk, &vrho, &vsigma,
           &v2rho2, &v2rhosigma, &v2sigma2, &v3rho3, &v3rho2sigma, &v3rhosigma2, &v3sigma3);

      }else{                                      /* polarized (general) case */
        func_pol(p, order, rho, sigma, OUT_PARAMS);
      } /* polarization */
    }
    
    internal_counters_gga_next
      (&(p->dim), 0, &rho, &sigma, &zk, &vrho, &vsigma,
       &v2rho2, &v2rhosigma, &v2sigma2, &v3rho3, &v3rho2sigma, &v3rhosigma2, &v3sigma3);
  }   /* for(ip) */
}
