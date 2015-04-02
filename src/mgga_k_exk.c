/*
 Copyright (C) 2015 D. Strubbe

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

#define XC_MGGA_K_EXK          666 /* exact kinetic energy */

static void work_mgga_k (const XC(func_type) *p, int np,
			 const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
			 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau,
			 FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
			 FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
			 FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau)
{
  XC(mgga_work_x_t) r;
  FLOAT dens;
  int ip;
   
  XC(rho2dzeta)(p->nspin, rho, &dens, &r.zeta);

  for(ip = 0; ip < np; ip++){
    if(zk          != NULL) { *zk          = *vtau; }
    if(vtau        != NULL) { *vtau        = 1.0;   }

    /* everything else is zero */
    if(vrho        != NULL) { *vrho        = 0.0;   }
    if(vsigma      != NULL) { *vsigma      = 0.0;   }
    if(vlapl       != NULL) { *vlapl       = 0.0;   }
    if(v2rho2      != NULL) { *v2rho2      = 0.0;   }
    if(v2sigma2    != NULL) { *v2sigma2    = 0.0;   }
    if(v2lapl2     != NULL) { *v2lapl2     = 0.0;   }
    if(v2tau2      != NULL) { *v2tau2      = 0.0;   }
    if(v2rhosigma  != NULL) { *v2rhosigma  = 0.0;   }
    if(v2rholapl   != NULL) { *v2rholapl   = 0.0;   }
    if(v2rhotau    != NULL) { *v2rhotau    = 0.0;   }
    if(v2sigmalapl != NULL) { *v2sigmalapl = 0.0;   }
    if(v2sigmatau  != NULL) { *v2sigmatau  = 0.0;   }
    if(v2lapltau   != NULL) { *v2lapltau   = 0.0;   }

    if(zk != NULL)
      *zk /= dens; /* we want energy per particle */

  end_ip_loop:
    /* increment pointers */
    rho   += p->n_rho;
    sigma += p->n_sigma;
    tau   += p->n_tau;
    lapl  += p->n_lapl;

    if(zk != NULL)
      zk += p->n_zk;

    if(vrho != NULL){
      vrho   += p->n_vrho;
      vsigma += p->n_vsigma;
      vtau   += p->n_vtau;
      vlapl  += p->n_vlapl;
    }

    if(v2rho2 != NULL){
      v2rho2      += p->n_v2rho2;
      v2sigma2    += p->n_v2sigma2;
      v2tau2      += p->n_v2tau2;
      v2lapl2     += p->n_v2lapl2;
      v2rhosigma  += p->n_v2rhosigma;
      v2rhotau    += p->n_v2rhotau;
      v2rholapl   += p->n_v2rholapl;
      v2sigmatau  += p->n_v2sigmatau;
      v2sigmalapl += p->n_v2sigmalapl;
      v2lapltau   += p->n_v2lapltau;
    }
  }
}

const XC(func_info_type) XC(func_info_mgga_k_exk) = {
  XC_MGGA_K_EXK,
  XC_KINETIC,
  "exact kinetic energy",
  XC_FAMILY_MGGA,
  {&xc_ref_Strubbe, &xc_ref_Schrodinger, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  0.0, 0.0, 0.0, 0.0,
  NULL,
  NULL, NULL,
  NULL,
  work_mgga_k
};
