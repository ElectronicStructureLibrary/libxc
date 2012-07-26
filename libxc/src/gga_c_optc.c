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

#define XC_GGA_C_OPTC       200 /* Optimized correlation functional of Cohen and Handy */

static void 
gga_c_optc_init(XC(func_type) *p)
{
  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  /* PW91 has always to be called polarized */
  XC(func_init)(p->func_aux[0], XC_GGA_C_PW91, XC_POLARIZED);
}

typedef struct XC(gga_work_c_t) {
  int   order; /* to which order should I return the derivatives */
  FLOAT rs, zeta, xt, xs[2];
  FLOAT f;
  FLOAT dfdrs, dfdz, dfdxt, dfdxs[2];
  FLOAT d2fdrs2, d2fdrsz, d2fdrsxt, d2fdrsxs[2], d2fdz2, 
    d2fdzxt, d2fdzxs[2], d2fdxt2, d2fdxtxs[2], d2fdxs2[3];
} XC(gga_work_c_t);

static inline void 
func(const XC(func_type) *p, int order, FLOAT rs, FLOAT zeta, FLOAT xt, FLOAT *xs,
     FLOAT *f, FLOAT *dfdrs, FLOAT *dfdz, FLOAT *dfdxt, FLOAT *dfdxs,
     FLOAT *d2fdrs2, FLOAT *d2fdrsz, FLOAT *d2fdrsxt, FLOAT *d2fdrsxs, FLOAT *d2fdz2, 
     FLOAT *d2fdzxt, FLOAT *d2fdzxs, FLOAT *d2fdxt2, FLOAT *d2fdxtxs, FLOAT *d2fdxs2)
{
  static FLOAT c1 = 1.1015, c2 = 0.6625;
  FLOAT opz, omz, copz, comz;

  XC(gga_work_c_t) f_par[2], f_anti;

  opz = 1.0 + zeta;
  omz = 1.0 - zeta;
  copz = CBRT(opz);
  comz = CBRT(omz);

  /* calculate the total part */
  f_anti.order = order;
  f_anti.rs    = rs;
  f_anti.zeta  = zeta;
  f_anti.xt    = xt;
  f_anti.xs[0] = xs[0];
  f_anti.xs[1] = xs[1];

  XC(gga_c_pw91_func) 
    (p->func_aux[0],  f_anti.order, f_anti.rs, f_anti.zeta, f_anti.xt, f_anti.xs,
     &(f_anti.f), &(f_anti.dfdrs), &(f_anti.dfdz), &(f_anti.dfdxt), f_anti.dfdxs,
     &(f_anti.d2fdrs2), &(f_anti.d2fdrsz), &(f_anti.d2fdrsxt), f_anti.d2fdrsxs, &(f_anti.d2fdz2), 
     (&f_anti.d2fdzxt), f_anti.d2fdzxs, &(f_anti.d2fdxt2), f_anti.d2fdxtxs, f_anti.d2fdxs2);

  /* now the spin up */
  f_par[0].order = order;
  f_par[0].rs    = rs*M_CBRT2*copz;
  f_par[0].zeta  = 1.0;
  f_par[0].xt    = xs[0];
  f_par[0].xs[0] = xs[0];
  f_par[0].xs[1] = 0.0;

  XC(gga_c_pw91_func) 
    (p->func_aux[0],  f_par[0].order, f_par[0].rs, f_par[0].zeta, f_par[0].xt, f_par[0].xs,
     &(f_par[0].f), &(f_par[0].dfdrs), &(f_par[0].dfdz), &(f_par[0].dfdxt), f_par[0].dfdxs,
     &(f_par[0].d2fdrs2), &(f_par[0].d2fdrsz), &(f_par[0].d2fdrsxt), f_par[0].d2fdrsxs, &(f_par[0].d2fdz2), 
     (&f_par[0].d2fdzxt), f_par[0].d2fdzxs, &(f_par[0].d2fdxt2), f_par[0].d2fdxtxs, f_par[0].d2fdxs2);

  /* now the spin down */
  f_par[1].order = order;
  f_par[1].rs    = rs*M_CBRT2*comz;
  f_par[1].zeta  = -1.0;
  f_par[1].xt    = xs[1];
  f_par[1].xs[0] = 0.0;
  f_par[1].xs[1] = xs[1];

  XC(gga_c_pw91_func) 
    (p->func_aux[0],  f_par[1].order, f_par[1].rs, f_par[1].zeta, f_par[1].xt, f_par[1].xs,
     &(f_par[1].f), &(f_par[1].dfdrs), &(f_par[1].dfdz), &(f_par[1].dfdxt), f_par[1].dfdxs,
     &(f_par[1].d2fdrs2), &(f_par[1].d2fdrsz), &(f_par[1].d2fdrsxt), f_par[1].d2fdrsxs, &(f_par[1].d2fdz2), 
     (&f_par[1].d2fdzxt), f_par[1].d2fdzxs, &(f_par[1].d2fdxt2), f_par[1].d2fdxtxs, f_par[1].d2fdxs2);

  /* now we add everything */
  
  *f = c1*f_anti.f + (c2 - c1)*(f_par[0].f + f_par[1].f);

  if(order < 1) return;

  *dfdrs   = c1*f_anti.dfdrs + (c2 - c1)*M_CBRT2*(f_par[0].dfdrs*copz + f_par[1].dfdrs*comz);
  *dfdz    = c1*f_anti.dfdz  + (c2 - c1)*M_CBRT2*rs/3.0*(f_par[0].dfdrs/(copz*copz) - f_par[1].dfdrs/(comz*comz));
  *dfdxt   = c1*f_anti.dfdxt;
  dfdxs[0] = c1*f_anti.dfdxs[0] + (c2 - c1)*(f_par[0].dfdxt + f_par[0].dfdxs[0]);
  dfdxs[1] = c1*f_anti.dfdxs[1] + (c2 - c1)*(f_par[1].dfdxt + f_par[1].dfdxs[1]);

  if(order < 2) return;

  *d2fdrs2    = c1*f_anti.d2fdrs2 + (c2 - c1)*M_CBRT2*M_CBRT2*(f_par[0].d2fdrs2*copz*copz + f_par[1].d2fdrs2*comz*comz);
  *d2fdrsz    = c1*f_anti.d2fdrsz + (c2 - c1)*M_CBRT2/3.0*
    (f_par[0].dfdrs/(copz*copz) - f_par[1].dfdrs/(comz*comz) +
     M_CBRT2*(f_par[0].d2fdrs2/copz - f_par[1].d2fdrs2/comz));
  *d2fdrsxt   = c1*f_anti.d2fdrsxt;
  d2fdrsxs[0] = c1*f_anti.d2fdrsxs[0] + (c2 - c1)*M_CBRT2*(f_par[0].d2fdrsxt + f_par[0].d2fdrsxs[0])*copz;
  d2fdrsxs[1] = c1*f_anti.d2fdrsxs[1] + (c2 - c1)*M_CBRT2*(f_par[1].d2fdrsxt + f_par[1].d2fdrsxs[1])*comz;
  *d2fdz2     = c1*f_anti.d2fdz2 + (c2 - c1)*M_CBRT2*rs/3.0*
    (-2.0/3.0*(f_par[0].dfdrs/(opz*copz*copz) + f_par[1].dfdrs/(omz*comz*comz)) + 
     M_CBRT2*rs/3.0*(f_par[0].d2fdrs2/(opz*copz) + f_par[1].d2fdrs2/(omz*comz)));
  *d2fdzxt    = c1*f_anti.d2fdzxt;
  d2fdzxs[0]  = c1*f_anti.d2fdzxs[0]  + (c2 - c1)*M_CBRT2*rs/3.0*(f_par[0].d2fdrsxt + f_par[0].d2fdrsxs[0])/(copz*copz);
  d2fdzxs[1]  = c1*f_anti.d2fdzxs[1]  + (c2 - c1)*M_CBRT2*rs/3.0*(f_par[1].d2fdrsxt + f_par[1].d2fdrsxs[1])/(comz*comz);
  *d2fdxt2    = c1*f_anti.d2fdxt2;
  d2fdxtxs[0] = c1*f_anti.d2fdxtxs[0];
  d2fdxtxs[1] = c1*f_anti.d2fdxtxs[1];
  d2fdxs2[0]  = c1*f_anti.d2fdxs2[0]  + (c2 - c1)*(f_par[0].d2fdxt2 + 2.0*f_par[0].d2fdxtxs[0] + f_par[0].d2fdxs2[0]);
  d2fdxs2[1]  = c1*f_anti.d2fdxs2[1];
  d2fdxs2[2]  = c1*f_anti.d2fdxs2[2]  + (c2 - c1)*(f_par[1].d2fdxt2 + 2.0*f_par[1].d2fdxtxs[1] + f_par[1].d2fdxs2[2]);
}

#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_c_optc) = {
  XC_GGA_C_OPTC,
  XC_CORRELATION,
  "Optimized correlation functional of Cohen and Handy",
  XC_FAMILY_GGA,
  "AJ Cohen and NC Handy, Mol. Phys. 99, 607-615 (2001)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-27, 1e-32, 0.0, 1e-32,
  gga_c_optc_init,
  NULL, NULL,
  work_gga_c,
};
