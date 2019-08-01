/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_mgga_deorb.c
 * @brief This file is to be included in deorbitalized MGGA functionals.
 */
   
#ifdef XC_NO_EXC
#define OUT_PARAMS(P_) MGGA_OUT_PARAMS_NO_EXC(P_)
#else
#define OUT_PARAMS(P_) P_ zk, MGGA_OUT_PARAMS_NO_EXC(P_)
#endif

/**
 * @param[in,out] func_type: pointer to functional structure
 */
static void 
work_mgga_deorb(const XC(func_type) *p, int np,
         const double *rho, const double *sigma, const double *lapl, const double *tau,
         OUT_PARAMS(double *))
{
  int ip, order;
  double dens, zeta;
  double *taus, *dtausdrho, *dtausdsigma, *dtausdlapl;
  const xc_dimensions *dim = &(p->dim);

  taus = (double *) malloc(sizeof(double)*np*dim->tau);
  dtausdrho = (double *) malloc(sizeof(double)*np*dim->vrho);
  dtausdsigma = (double *) malloc(sizeof(double)*np*dim->vsigma);
  dtausdlapl = (double *) malloc(sizeof(double)*np*dim->vlapl);

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
        ked_unpol(order, rho, sigma, lapl, taus, dtausdrho, dtausdsigma, dtausdlapl);
        func_unpol(p, order, rho, sigma, lapl, taus, OUT_PARAMS());

        /* Deorbitalization Contribution */
        vrho[0] = vrho[0] + vtau[0] * dtausdrho[0];
        vsigma[0] = vsigma[0] + vtau[0] * dtausdsigma[0];
        vlapl[0] = vlapl[0] + vtau[0] * dtausdlapl[0];
        vtau[0] = 0.0;
      
      }else if(zeta >  1.0 - 1e-10){              /* ferromagnetic case - spin 0 */
        ked_ferr(order, rho, sigma, lapl, taus, dtausdrho, dtausdsigma, dtausdlapl);
        func_ferr(p, order, rho, sigma, lapl, taus, OUT_PARAMS());
        /* Deorbitalization contribution */
        vrho[0] = vrho[0] + vtau[0] * dtausdrho[0];
        vsigma[0] = vsigma[0] + vtau[0] * dtausdsigma[0];
        vlapl[0] = vlapl[0] + vtau[0] * dtausdlapl[0];
        vtau[0] = 0.0;

        
      }else if(zeta < -1.0 + 1e-10){              /* ferromagnetic case - spin 1 */
        internal_counters_mgga_next(&(p->dim), -1, &rho, &sigma, &lapl, &tau, &zk, MGGA_OUT_PARAMS_NO_EXC(&));
        ked_ferr(order, rho, sigma, lapl, taus, dtausdrho, dtausdsigma, dtausdlapl);
        func_ferr(p, order, rho, sigma, lapl, taus, OUT_PARAMS());
        /* Deorbitalization contribution */
        vrho[0] = vrho[0] + vtau[0] * dtausdrho[0];
        vsigma[0] = vsigma[0] + vtau[0] * dtausdsigma[0];
        vlapl[0] = vlapl[0] + vtau[0] * dtausdlapl[0];
        vtau[0] = 0.0;
        internal_counters_mgga_prev(&(p->dim), -1, &rho, &sigma, &lapl, &tau, &zk, MGGA_OUT_PARAMS_NO_EXC(&));
      }else{                                      /* polarized (general) case */
        ked_pol(order, rho, sigma, lapl, taus, dtausdrho, dtausdsigma, dtausdlapl);
        func_pol(p, order, rho, sigma, lapl, taus, OUT_PARAMS());
        /* Deorbitalization contribution */
        vrho[0] = vrho[0] + vtau[0] * dtausdrho[0];
        vsigma[0] = vsigma[0] + vtau[0] * dtausdsigma[0];
        vlapl[0] = vlapl[0] + vtau[0] * dtausdlapl[0];
        vtau[0] = 0.0;

        vrho[1] = vrho[1] + vtau[1] * dtausdrho[1];
        vsigma[2] = vsigma[2] + vtau[1] * dtausdsigma[2];
        vlapl[1] = vlapl[1] + vtau[1] * dtausdlapl[1];
        vtau[1] = 0.0;

      } /* polarization */
    }
    
    internal_counters_mgga_next(&(p->dim), 0, &rho, &sigma, &lapl, &tau, &zk, MGGA_OUT_PARAMS_NO_EXC(&));
  }   /* for(ip) */
};
