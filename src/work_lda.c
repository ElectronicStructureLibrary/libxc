/*
 Copyright (C) 2006-2018 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_lda.c
 * @brief This file is to be included in LDA functionals.
 */


/* hack to avoid compiler warnings */
#define NOARG

#ifdef XC_NO_EXC
#define OUT_PARAMS LDA_OUT_PARAMS_NO_EXC(NOARG)
#else
#define OUT_PARAMS zk, LDA_OUT_PARAMS_NO_EXC(NOARG)
#endif

/**
 * @param[in,out] func_type: pointer to functional structure
 */
static void 
work_lda(const XC(func_type) *p, size_t np, const double *rho, 
	 double *zk, LDA_OUT_PARAMS_NO_EXC(double *))
{
  size_t ip;
  int order;
  double dens, zeta;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  for(ip = 0; ip < np; ip++){
    xc_rho2dzeta(p->nspin, rho, &dens, &zeta);

    if(dens > p->dens_threshold){
      if(p->nspin == XC_UNPOLARIZED){             /* unpolarized case */
        func_unpol(p, order, rho, OUT_PARAMS);
      
      }else if(zeta >  1.0 - 1e-10){              /* ferromagnetic case - spin 0 */
        func_ferr(p, order, rho, OUT_PARAMS);
        
      }else if(zeta < -1.0 + 1e-10){              /* ferromagnetic case - spin 1 */
        internal_counters_lda_next(&(p->dim), -1, &rho, &zk, LDA_OUT_PARAMS_NO_EXC(&));
        func_ferr(p, order, rho, OUT_PARAMS);
        internal_counters_lda_prev(&(p->dim), -1, &rho, &zk, LDA_OUT_PARAMS_NO_EXC(&));

      }else{                                      /* polarized (general) case */
        func_pol(p, order, rho, OUT_PARAMS);
      } /* polarization */
    }
    
    internal_counters_lda_next(&(p->dim), 0, &rho, &zk, LDA_OUT_PARAMS_NO_EXC(&));
  }   /* for(ip) */
}
