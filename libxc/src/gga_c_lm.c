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

/************************************************************************
  This functional is provided for historical reasons.
  It was one of the first GGAs that ever appeared.
************************************************************************/

#define XC_GGA_C_LM          137 /* Langreth and Mehl correlation          */

static void 
gga_c_lm_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_vBH, p->nspin);
}


static inline void 
func(const XC(gga_type) *p, int order, FLOAT rs, FLOAT zeta, FLOAT xt, FLOAT *xs,
     FLOAT *f, FLOAT *dfdrs, FLOAT *dfdz, FLOAT *dfdxt, FLOAT *dfdxs,
     FLOAT *d2fdrs2, FLOAT *d2fdrsz, FLOAT *d2fdrsxt, FLOAT *d2fdrsxs, FLOAT *d2fdz2, 
     FLOAT *d2fdzxt, FLOAT *d2fdzxs, FLOAT *d2fdxt2, FLOAT *d2fdxtxs, FLOAT *d2fdxs2)
{
  const FLOAT a1 = 4.28e-3/2.0; /* The factor of 2 converts from Rydberg to Hartree */
  const FLOAT a2 = -0.262;
  const FLOAT a3 = 7.0/9.0;

  FLOAT alpha;
  FLOAT H, dHdx1, d2Hdx12, dx1drs, dx1dxt, d2x1drs2, d2x1drsxt;
  FLOAT x1;

  XC(lda_rs_zeta) pw;

  alpha = POW(4.0*M_PI/3.0, 1.0/6.0);

  pw.order = order;
  pw.rs[0] = SQRT(rs);
  pw.rs[1] = rs;
  pw.rs[2] = rs*rs;
  pw.zeta  = zeta;

  XC(lda_c_hl_func)(p->func_aux[0]->lda, &pw);

  x1 = xt/(alpha*pw.rs[0]);
  
  H = a1*x1*x1*(exp(a2*x1) - a3);
  *f = pw.zk + H;

  if(order < 1) return;

  dx1drs = -xt/(2.0*alpha*rs*pw.rs[0]);
  dx1dxt = 1.0/(alpha*pw.rs[0]);
  dHdx1  = a1*x1*((2.0 + a2*x1)*exp(a2*x1) - 2.0*a3);

  *dfdrs   = pw.dedrs + dHdx1*dx1drs;
  *dfdz    = pw.dedz;
  *dfdxt   = dHdx1*dx1dxt;
  dfdxs[0] = 0.0;
  dfdxs[1] = 0.0;

  if(order < 2) return;

  d2x1drs2  = 3.0*xt/(4.0*alpha*pw.rs[2]*pw.rs[0]);
  d2x1drsxt = -1.0/(2.0*alpha*rs*pw.rs[0]);
  d2Hdx12   = a1*((2.0 + a2*x1*(4.0 + a2*x1))*exp(a2*x1) - 2.0*a3);

  *d2fdrs2    = pw.d2edrs2 + d2Hdx12*dx1drs*dx1drs + dHdx1*d2x1drs2;
  *d2fdrsz    = pw.d2edrsz;
  *d2fdrsxt   = d2Hdx12*dx1drs*dx1dxt + dHdx1*d2x1drsxt;
  d2fdrsxs[0] = 0.0;
  d2fdrsxs[1] = 0.0;
  *d2fdz2     = pw.d2edz2;
  *d2fdzxt    = 0.0;
  d2fdzxs[0]  = 0.0;
  d2fdzxs[1]  = 0.0;
  *d2fdxt2    = d2Hdx12*dx1dxt*dx1dxt;
  d2fdxtxs[0] = 0.0;
  d2fdxtxs[1] = 0.0;
  d2fdxs2[0]  = 0.0;
  d2fdxs2[1]  = 0.0;
  d2fdxs2[2]  = 0.0;
}

#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_lm) = {
  XC_GGA_C_LM,
  XC_CORRELATION,
  "Langreth & Mehl",
  XC_FAMILY_GGA,
  "DC Langreth and MJ Mehl, Phys. Rev. Lett. 47, 446 (1981)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_c_lm_init,
  NULL,
  NULL,            /* this is not an LDA                   */
  work_gga_c,
};

