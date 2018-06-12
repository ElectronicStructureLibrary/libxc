/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


static void 
work_mgga_c
(const xc_func_type *p, int np,
 const double *rho, const double *sigma, const double *lapl, const double *tau,
 double *zk,
 double *vrho, double *vsigma, double *vlapl, double *vtau,
 double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
 double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
 double *v2lapl2,  double *v2taulapl,
 double *v2tau2,
 double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
 double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
 double *v3rholapl2, double *v3rhotaulapl,
 double *v3rhotau2,
 double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
 double *v3sigmalapl2, double *v3sigmataulapl,
 double *v3sigmatau2,
 double *v3lapl3,  double *v3taulapl2,
 double *v3tau2lapl,
 double *v3tau3
)
{
  xc_mgga_work_c_t r;
  double min_grad2 = p->dens_threshold*p->dens_threshold, min_tau = p->dens_threshold;
  int ip;

  /* set all elements of r to zero */
  memset(&r, 0, sizeof(r));

  r.order = -1;
  if(zk     != NULL) r.order = 0;
  if(vrho   != NULL) r.order = 1;
  if(v2rho2 != NULL) r.order = 2;

  if(r.order < 0) return;

  for(ip = 0; ip < np; ip++){
    double rho13[3], drs, dxt;
    double ndzdn[2], dxsdn[2];
    double dxtds, dxsds[2];
    double dusdn[2], dusdlapl[2], dtsdn[2], dtsdtau[2];

    xc_rho2dzeta(p->nspin, rho, &(r.dens), &(r.z));

    if(r.dens < p->dens_threshold) goto end_ip_loop;
    
    r.rs = RS(r.dens);
    rho13[2] = CBRT(r.dens);

    if(p->nspin == XC_UNPOLARIZED){
      r.ds[0]  = r.dens/2.0;
      r.ds[1]  = r.ds[0];

      rho13[0] = rho13[2]/M_CBRT2;
      rho13[1] = rho13[0];

      /* we already know that dens > min_dens */
      r.sigmat = max(min_grad2, sigma[0]);
      r.xt     = sqrt(r.sigmat)/(r.dens*rho13[2]);

      r.sigmas[0] = r.sigmat/4.0;
      r.sigmas[1] = r.sigmas[0];
      r.sigmas[2] = r.sigmas[0];

      r.xs[0]  = M_CBRT2*r.xt;
      r.xs[1]  = r.xs[0];

      if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
        r.us[0]  = lapl[0]/(2.0*r.ds[0]*rho13[0]*rho13[0]); /* lapl/rho^(5/3) */
        r.us[1]  = r.us[0];
      }

      r.ts[0]  = max(2.0*min_tau, tau[0]/(2.0*r.ds[0]*rho13[0]*rho13[0]));  /* tau/rho^(5/3) */
      r.ts[1]  = r.ts[0];
    }else{
      /* there are lots of derivatives that involve inverse
         powers of (1 +- z). For these not to give NaN, we
         must have abs(1 +- z) > DBL_EPSILON                 */
      if(1.0 + r.z < DBL_EPSILON) r.z = -1.0 + DBL_EPSILON;
      if(1.0 - r.z < DBL_EPSILON) r.z =  1.0 - DBL_EPSILON;

      r.ds[0]  = max(p->dens_threshold, rho[0]);
      r.ds[1]  = max(p->dens_threshold, rho[1]);

      rho13[0] = CBRT(r.ds[0]);
      rho13[1] = CBRT(r.ds[1]);
      
      r.sigmat = max(min_grad2, sigma[0] + 2.0*sigma[1] + sigma[2]);
      r.xt     = sqrt(r.sigmat)/(r.dens*rho13[2]);
      
      r.sigmas[0] = max(min_grad2, sigma[0]);
      r.sigmas[1] = max(min_grad2, sigma[1]);
      r.sigmas[2] = max(min_grad2, sigma[2]);

      r.xs[0] = sqrt(r.sigmas[0])/(r.ds[0]*rho13[0]);
      r.xs[1] = sqrt(r.sigmas[2])/(r.ds[1]*rho13[1]);

      if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
        r.us[0]   = lapl[0]/(r.ds[0]*rho13[0]*rho13[0]);
        r.us[1]   = lapl[1]/(r.ds[1]*rho13[1]*rho13[1]);
      }

      r.ts[0]   = max(min_tau, tau[0]/(r.ds[0]*rho13[0]*rho13[0]));
      r.ts[1]   = max(min_tau, tau[1]/(r.ds[1]*rho13[1]*rho13[1]));
    }
  
    func(p, &r);

    if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
      *zk = r.f;

    if(r.order < 1) goto end_ip_loop;
    
    /* setup auxiliary variables */
    drs   =     -r.rs/(3.0*r.dens);
    dxt   = -4.0*r.xt/(3.0*r.dens);
    dxtds = r.xt/(2.0*r.sigmat);

    if(p->nspin == XC_POLARIZED){
      ndzdn[1]    = -(r.z + 1.0);
      ndzdn[0]    = -(r.z - 1.0);

      dxsdn[1]    = -4.0*r.xs[1]/(3.0*r.ds[1]);
      dxsdn[0]    = -4.0*r.xs[0]/(3.0*r.ds[0]);

      dxsds[1]    = r.xs[1]/(2.0*r.sigmas[2]);
      dxsds[0]    = r.xs[0]/(2.0*r.sigmas[0]);

      dtsdn[1]    = -5.0*r.ts[1]/(3.0*r.ds[1]);
      dtsdn[0]    = -5.0*r.ts[0]/(3.0*r.ds[0]);

      dtsdtau[1]  = 1.0/(r.ds[1]*rho13[1]*rho13[1]);
      dtsdtau[0]  = 1.0/(r.ds[0]*rho13[0]*rho13[0]);

      if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
        dusdn[1]    = -5.0*r.us[1]/(3.0*r.ds[1]);
        dusdn[0]    = -5.0*r.us[0]/(3.0*r.ds[0]);

        dusdlapl[1] = dtsdtau[1];
        dusdlapl[0] = dtsdtau[0];
      }

    }else{
      dxsdn[0]    = M_CBRT2*dxt;
      dxsds[0]    = M_CBRT2*dxtds;

      dtsdn[0]    = -5.0*r.ts[0]/(6.0*r.ds[0]);
      dtsdtau[0]  = 1.0/(2.0*r.ds[0]*rho13[0]*rho13[0]);

      if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
        dusdn[0]    = -5.0*r.us[0]/(6.0*r.ds[0]);
        dusdlapl[0] = dtsdtau[0];
      } 
    }

    if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
      vrho[0]   = r.f + r.dens*(r.dfdrs*drs + r.dfdxt*dxt);
      vsigma[0] = r.dens*r.dfdxt*dxtds;

      if(p->nspin == XC_POLARIZED){
	vrho[1]   = vrho[0] + r.dfdz*ndzdn[1] + r.dens*(r.dfdxs[1]*dxsdn[1] + r.dfdts[1]*dtsdn[1]);
	vrho[0]   = vrho[0] + r.dfdz*ndzdn[0] + r.dens*(r.dfdxs[0]*dxsdn[0] + r.dfdts[0]*dtsdn[0]);

	vsigma[2] = vsigma[0] + r.dens*r.dfdxs[1]*dxsds[1];
	vsigma[1] = 2.0*vsigma[0];
	vsigma[0] = vsigma[0] + r.dens*r.dfdxs[0]*dxsds[0];

        if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
          vrho[1]  += r.dens*r.dfdus[1]*dusdn[1];
          vrho[0]  += r.dens*r.dfdus[0]*dusdn[0];

          vlapl[1]  = r.dens*r.dfdus[1]*dusdlapl[1];
          vlapl[0]  = r.dens*r.dfdus[0]*dusdlapl[0];
        }

	vtau[1]   = r.dens*r.dfdts[1]*dtsdtau[1];
	vtau[0]   = r.dens*r.dfdts[0]*dtsdtau[0];
	
      }else{
	 /* factor of 2 comes from sum over sigma */
	vrho[0]   += 2.0*r.dens*(r.dfdxs[0]*dxsdn[0] + r.dfdts[0]*dtsdn[0]);
	vsigma[0] += 2.0*r.dens*r.dfdxs[0]*dxsds[0];

        if(p->info->flags & XC_FLAGS_NEEDS_LAPLACIAN){
          vrho[0]   += 2.0*r.dens*r.dfdus[0]*dusdn[0];

          vlapl[0]   = 2.0*r.dens*r.dfdus[0]*dusdlapl[0];
        }

	vtau[0]    = 2.0*r.dens*r.dfdts[0]*dtsdtau[0];
      }
    }
    
    if(r.order < 2) goto end_ip_loop;

  end_ip_loop:
internal_counters_mgga_next(&(p->dim), 0, &rho, &sigma, &lapl, &tau,
                                &zk, &vrho, &vsigma, &vlapl, &vtau,
                                &v2rho2, &v2rhosigma, &v2rholapl, &v2rhotau,
                                &v2sigma2, &v2sigmalapl, &v2sigmatau,
                                &v2lapl2, &v2taulapl,
                                &v2tau2,
                                &v3rho3, &v3rho2sigma, &v3rho2lapl, &v3rho2tau,
                                &v3rhosigma2, &v3rhosigmalapl, &v3rhosigmatau,
                                &v3rholapl2, &v3rhotaulapl,
                                &v3rhotau2,
                                &v3sigma3, &v3sigma2lapl, &v3sigma2tau,
                                &v3sigmalapl2, &v3sigmataulapl,
                                &v3sigmatau2,
                                &v3lapl3, &v3taulapl2,
                                &v3tau2lapl,
                                &v3tau3);
  }
}
