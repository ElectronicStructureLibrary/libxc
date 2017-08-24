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

#include "util.h"

#define XC_GGA_X_ITYH 529 /* short-range recipe for exchange GGA functionals */

typedef struct{
  int func_id;
  xc_gga_enhancement_t enhancement_factor;
} gga_x_ityh_params;

static void
gga_x_ityh_init(XC(func_type) *p)
{
  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_ityh_params));

  /* random functional, mainly intended for testing */
  ((gga_x_ityh_params *) (p->params))->func_id = -1;
  XC(gga_x_ityh_set_params)(p, XC_GGA_X_B88, 0.2);
}

void 
XC(gga_x_ityh_set_params)(XC(func_type) *p, int func_id, FLOAT omega)
{
  gga_x_ityh_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ityh_params *) (p->params);

  p->cam_omega = omega;

  /* if func_id == -1 do nothing */
  if(func_id != -1 && params->func_id == -1){ /* intialize stuff */
    p->n_func_aux  = 1;
    p->func_aux    = (XC(func_type) **) malloc(sizeof(XC(func_type) *));
    p->func_aux[0] = (XC(func_type)  *) malloc(sizeof(XC(func_type)  ));
  }

  if(func_id != -1 && params->func_id != func_id){
    if(params->func_id != -1)
      XC(func_end) (p->func_aux[0]);

    params->func_id = func_id;
    XC(func_init) (p->func_aux[0], params->func_id, p->nspin);

    params->enhancement_factor = XC(get_gga_enhancement_factor)(func_id);
  }
}


static inline void 
func_3(const XC(func_type) *p, int order, FLOAT x, FLOAT ds,
     FLOAT *f, FLOAT *dfdx, FLOAT *lvrho)
{
  gga_x_ityh_params *params;
  XC(gga_work_x_t) aux;

  FLOAT k_GGA, K_GGA, aa, f_aa, df_aa, d2f_aa, d3f_aa;
  FLOAT dk_GGAdr, dk_GGAdx, daadr, daadx;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_ityh_params *) (p->params);

  /* call enhancement factor */
  aux.x    = x;
  aux.order = order;
  params->enhancement_factor(p->func_aux[0], &aux);

  K_GGA = 2.0*X_FACTOR_C*aux.f;
  k_GGA = SQRT(9.0*M_PI/K_GGA)*CBRT(ds);

  aa = p->cam_omega/(2.0*k_GGA);

  XC(lda_x_attenuation_function)(XC_RSF_ERF, order, aa, &f_aa, &df_aa, &d2f_aa, &d3f_aa);

  *f = aux.f*f_aa;

  if(order < 1) return;

  dk_GGAdr =  k_GGA/(3.0*ds);
  dk_GGAdx = -k_GGA*aux.dfdx/(2.0*aux.f);

  daadr   = -aa*dk_GGAdr/k_GGA;
  daadx   = -aa*dk_GGAdx/k_GGA;

  *dfdx   = aux.dfdx*f_aa + aux.f*df_aa*daadx;
  *lvrho  = aux.f*df_aa*daadr; 
}

/* convert into work_gga_c_ variables */
static inline void 
func(const XC(func_type) *p, XC(gga_work_c_t) *r)
{
  int i;
  FLOAT ds, ex, f, lvrho, dexdrs, ddsdrs, dexdz, ddsdz;
  FLOAT sign[2] = {1.0, -1.0};

  r->f     = 0.0;
  r->dfdrs = 0.0;
  r->dfdz  = 0.0;

  for(i=0; i<2; i++){
    ds = POW(RS_FACTOR/r->rs, 3.0)*(1.0 + sign[i]*r->z)/2.0;
    func_3(p, r->order, r->xs[i], ds, &f, &(r->dfdxs[i]), &lvrho);

    ex = -X_FACTOR_C*RS_FACTOR*POW((1.0 + sign[i]*r->z)/2.0, 4.0/3.0)/r->rs;

    r->f += ex*f;

    if(r->order < 1) continue;

    ddsdrs = -3.0*ds/r->rs;
    dexdrs = -ex/r->rs;

    ddsdz  = POW(RS_FACTOR/r->rs, 3.0)*sign[i]*r->z/2.0;
    dexdz  = -4.0/6.0*sign[i]*X_FACTOR_C*RS_FACTOR*POW((1.0 + sign[i]*r->z)/2.0, 1.0/3.0)/r->rs;

    r->dfdrs    += dexdrs*f + ex*lvrho*ddsdrs;
    r->dfdz     += dexdz *f + ex*lvrho*ddsdz;
    r->dfdxs[i] *= ex;
  }
}

#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_x_ityh) = {
  XC_GGA_X_ITYH,
  XC_EXCHANGE,
  "Short-range recipe for exchange GGA functionals",
  XC_FAMILY_GGA,
  {&xc_ref_Iikura2001_3540, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  gga_x_ityh_init, NULL, 
  NULL, work_gga_c, NULL
};
