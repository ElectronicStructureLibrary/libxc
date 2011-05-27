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
#include <math.h>
#include "util.h"

static void 
work_gga_c(const void *p_, int np, const FLOAT *rho, const FLOAT *sigma,
	   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const XC(gga_type) *p = p_;

  FLOAT min_grad2 = MIN_GRAD*MIN_GRAD;
  int ip, order;

  order = -1;
  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  if(order < 0) return;

  for(ip = 0; ip < np; ip++){
    FLOAT dens, zeta, rs, ds[2], sigmat, sigmas[3], xt, xs[2];
    FLOAT f, drs, dxt, dfdrs, dfdz, dfdxt, dfdxs[2];
    FLOAT d2rs, d2xt, d2fdrs2, d2fdrsz, d2fdrsxt, d2fdrsxs[2], d2fdz2, d2fdzxt, d2fdzxs[2], d2fdxt2, d2fdxtxs[2], d2fdxs2[3];
    FLOAT ndzdn[2], dxsdn[2];
    FLOAT dxtds, d2xtds2, d2xtdns, dxsds[2], d2xsds2[2], d2xsdns[2];

    XC(rho2dzeta)(p->nspin, rho, &dens, &zeta);

    if(dens < MIN_DENS) goto end_ip_loop;

    rs = RS(dens);
    if(p->nspin == XC_UNPOLARIZED){
      ds[0]  = dens/2.0;
      ds[1]  = ds[0];

      sigmat = max(min_grad2, sigma[0]);
      xt     = SQRT(sigmat)/ POW(dens, 4.0/3.0);

      sigmas[0] = sigmat/4.0;
      sigmas[1] = sigmas[0];
      sigmas[2] = sigmas[0];

      xs[0]  = CBRT(2.0)*xt;
      xs[1]  = xs[0];
    }else{
      ds[0]  = max(MIN_DENS, rho[0]);
      ds[1]  = max(MIN_DENS, rho[1]);
      
      sigmat = max(min_grad2, sigma[0] + 2.0*sigma[1] + sigma[2]);
      xt     = SQRT(sigmat)/ POW(dens, 4.0/3.0);
      
      sigmas[0] = max(min_grad2, sigma[0]);
      sigmas[1] = max(min_grad2, sigma[1]);
      sigmas[2] = max(min_grad2, sigma[2]);

      xs[0] = SQRT(sigmas[0])/POW(ds[0], 4.0/3.0);
      xs[1] = SQRT(sigmas[2])/POW(ds[1], 4.0/3.0);
    }
  
    func(p, order, rs, zeta, xt, xs, &f, &dfdrs, &dfdz, &dfdxt, dfdxs,
	 &d2fdrs2, &d2fdrsz, &d2fdrsxt, d2fdrsxs, &d2fdz2, &d2fdzxt, d2fdzxs, &d2fdxt2, d2fdxtxs, d2fdxs2);

    if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
      *zk = f;

    if(order < 1) goto end_ip_loop;
    
    drs =     -rs/(3.0*dens);
    dxt = -4.0*xt/(3.0*dens);

    if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC)){
      vrho[0]   = f + dens*(dfdrs*drs + dfdxt*dxt);

      dxtds = xt/(2.0*sigmat);
      vsigma[0] = dens*dfdxt*dxtds;

      if(p->nspin == XC_POLARIZED){
	/* vrho */
	ndzdn[1] = -(zeta + 1.0);
	ndzdn[0] = -(zeta - 1.0);

	dxsdn[1] = -4.0/3.0*xs[1]/ds[1];
	dxsdn[0] = -4.0/3.0*xs[0]/ds[0];

        vrho[1] = vrho[0] + dfdz*ndzdn[1] + dens*dfdxs[1]*dxsdn[1];
        vrho[0] = vrho[0] + dfdz*ndzdn[0] + dens*dfdxs[0]*dxsdn[0];;

	/* vsigma */
	dxsds[1] = xs[1]/(2.0*sigmas[2]);
	dxsds[0] = xs[0]/(2.0*sigmas[0]);

	vsigma[2] = vsigma[0] + dens*dfdxs[1]*dxsds[1];
	vsigma[1] = 2.0*vsigma[0];
	vsigma[0] = vsigma[0] + dens*dfdxs[0]*dxsds[0];

      }else{
	vrho[0]   += -8.0/3.0*xs[0]*dfdxs[0];
	vsigma[0] += dens*dfdxs[0]*xs[0]/(4.0*sigmas[0]);
      }
    }

    if(order < 2) goto end_ip_loop;
  
    d2rs = -4.0*drs/(3.0*dens);
    d2xt = -7.0*dxt/(3.0*dens);

    if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC)){
      v2rho2[0] = 2.0*dfdrs*drs + 2.0*dfdxt*dxt +
	dens*(d2fdrs2*drs*drs + d2fdxt2*dxt*dxt + dfdrs*d2rs + dfdxt*d2xt + 2.0*d2fdrsxt*drs*dxt);

      d2xtds2       = -dxtds/(2.0*sigmat);
      v2sigma2[0]   = dens*(d2fdxt2*dxtds*dxtds + dfdxt*d2xtds2);

      d2xtdns       = dxt/(2.0*sigmat);
      v2rhosigma[0] = dfdxt*dxtds + dens*(d2fdrsxt*drs*dxtds + d2fdxt2*dxt*dxtds + dfdxt*d2xtdns);

      if(p->nspin == XC_POLARIZED){
        FLOAT sign[3][2] = {{1.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}};
	int is;

        for(is=2; is>=0; is--){
	  int s1 = (is > 1) ?  1 :  0;         /* {0, 0, 1}[is] */
	  int s2 = (is > 0) ?  1 :  0;         /* {0, 1, 1}[is] */
	  FLOAT delta = (is == 1) ? 0.0 : 1.0; /* delta_{s1,s2} */

	  /* I have to simplify this! */

	  /* derivative wrt rs and xt */
	  v2rho2[is] = v2rho2[0] +
	    (- (zeta - sign[is][0])*(d2fdrsz*drs + d2fdzxt*dxt) 
	     -  4.0/3.0*xs[s1]*2.0/(1.0 + sign[is][0]*zeta)*(d2fdrsxs[s1]*drs + d2fdxtxs[s1]*dxt));

	  /* derivative wrt zeta */
          v2rho2[is] += -(zeta - sign[is][1])/dens *
	    (+ dens*(d2fdrsz*drs + d2fdzxt*dxt)
	     - (zeta - sign[is][0])*d2fdz2
	     - 4.0/3.0*xs[s1]*(d2fdzxs[s1] - dfdxs[s1]*sign[is][0]/(1.0 + sign[is][0]*zeta))*2.0/(1.0 + sign[is][0]*zeta)
	     );

	  /* derivative wrt xs */
	  v2rho2[is] += -4.0/3.0*xs[s2]*2.0/(1.0 + sign[is][1]*zeta)/dens *
	    (dfdxs[s2] + dens*(d2fdrsxs[s2]*drs + d2fdxtxs[s2]*dxt)
	     - (zeta - sign[is][0])*d2fdzxs[s2]
	     - 4.0/3.0*2.0/(1.0 + sign[is][0]*zeta)*(delta*dfdxs[s1] + xs[s1]*d2fdxs2[is])
	     );
	}

	/* v2sigma */
	d2xsds2[1] = -dxsds[1]/(2.0*sigmas[2]);
	d2xsds2[0] = -dxsds[0]/(2.0*sigmas[0]);

	v2sigma2[5] =     v2sigma2[0] + dens*
	  (2.0*d2fdxtxs[1]*dxtds*dxsds[1] + d2fdxs2[2]*dxsds[1]*dxsds[1] + dfdxs[1]*d2xsds2[1]);
	v2sigma2[4] = 2.0*v2sigma2[0] + dens*
	  (2.0*d2fdxtxs[1]*dxtds*dxsds[1]);
	v2sigma2[3] = 4.0*v2sigma2[0] + dens*
	  (4.0*d2fdxtxs[1]*dxtds*dxsds[1]);
	v2sigma2[2] =     v2sigma2[0] + dens*
	  (    d2fdxtxs[0]*dxtds*dxsds[2] + dxtds*(d2fdxtxs[0]*dxsds[1] + d2fdxtxs[1]*dxsds[0]) + d2fdxs2[1]*dxsds[0]*dxsds[1]);
	v2sigma2[1] =  2.0*v2sigma2[0] + dens*
	  (    d2fdxtxs[0]*dxtds*dxsds[1]);
	v2sigma2[0] =     v2sigma2[0] + dens*
	  (2.0*d2fdxtxs[0]*dxtds*dxsds[0] + d2fdxs2[0]*dxsds[0]*dxsds[0] + dfdxs[0]*d2xsds2[0]);

	/* v2rhosigma */
	d2xsdns[1] = -4.0/3.0*dxsds[1]/ds[1];
	d2xsdns[0] = -4.0/3.0*dxsds[0]/ds[0];

	v2rhosigma[5] =    v2rhosigma[0] + dfdxs[1]*dxsds[1] + ndzdn[1]*(d2fdzxt*dxtds + d2fdzxs[1]*dxsds[1]) + 
	  dens*(d2fdrsxs[1]*drs*dxsds[1] + d2fdxtxs[1]*dxsdn[1]*dxtds + d2fdxs2[2]*dxsdn[1]*dxsds[1] + dfdxs[1]*d2xsdns[1]);

	v2rhosigma[4] = 2.0*v2rhosigma[0] + 2.0*ndzdn[1]*d2fdzxt*dxtds + 
	  2.0*dens*d2fdxtxs[1]*dxsdn[1]*dxtds;

	v2rhosigma[3] =    v2rhosigma[0] + dfdxs[0]*dxsds[0] + ndzdn[1]*(d2fdzxt*dxtds + d2fdzxs[0]*dxsds[0]) + 
	  dens*(d2fdrsxs[0]*drs*dxsds[0] + d2fdxtxs[1]*dxsdn[1]*dxtds + d2fdxs2[1]*dxsdn[1]*dxsds[0]);

	v2rhosigma[2] =    v2rhosigma[0] + dfdxs[1]*dxsds[1] + ndzdn[0]*(d2fdzxt*dxtds + d2fdzxs[1]*dxsds[1]) + 
	  dens*(d2fdrsxs[1]*drs*dxsds[1] + d2fdxtxs[0]*dxsdn[0]*dxtds + d2fdxs2[1]*dxsdn[0]*dxsds[1]);

	v2rhosigma[1] = 2.0*v2rhosigma[0] + 2.0*ndzdn[0]*d2fdzxt*dxtds + 
	  2.0*dens*d2fdxtxs[0]*dxsdn[0]*dxtds;
	
	v2rhosigma[0] =    v2rhosigma[0] + dfdxs[0]*dxsds[0] + ndzdn[0]*(d2fdzxt*dxtds + d2fdzxs[0]*dxsds[0]) + 
	  dens*(d2fdrsxs[0]*drs*dxsds[0] + d2fdxtxs[0]*dxsdn[0]*dxtds + d2fdxs2[0]*dxsdn[0]*dxsds[0] + dfdxs[0]*d2xsdns[0]);

      }else{
	v2rho2[0] += -8.0*xs[0]/(3.0*dens) * (dfdxs[0] + 2.0*dens*(d2fdrsxs[0]*drs + d2fdxtxs[0]*dxt))
	  + (8.0*4.0)/(3.0*3.0)*(dfdxs[0] + 2.0*xs[0]*d2fdxs2[0])*xs[0]/dens;

	v2sigma2[0]   += dens*(2.0*d2fdxtxs[0]*xt*xs[0]/(8.0*sigmat*sigmas[0]) +
			       xs[0]/(32.0*sigmas[0]*sigmas[0])*(d2fdxs2[0]*xs[0] - dfdxs[0]));
	v2rhosigma[0] += xs[0]/(4.0*sigmas[0]) * 
	  (-4.0/3.0*(xt*d2fdxtxs[0] + dfdxs[0]/4.0 + xs[0]*d2fdxs2[0]) +
	   dens*(d2fdrsxs[0]*drs + d2fdxtxs[0]*dxt));
      }
    }

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
  }    
}
