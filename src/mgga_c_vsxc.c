/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_C_VSXC          232 /* VSxc from Van Voorhis and Scuseria (correlation part) */

typedef struct{
  const FLOAT alpha_ss, alpha_ab;
  const FLOAT dss[6], dab[6];
} mgga_c_vsxc_params;

static const mgga_c_vsxc_params par_vsxc = {
  0.00515088, 0.00304966,
  { 3.270912e-01, -3.228915e-02, -2.942406e-02,  2.134222e-03, -5.451559e-03,  1.577575e-02},
  { 7.035010e-01,  7.694574e-03,  5.152765e-02,  3.394308e-05, -1.269420e-03,  1.296118e-03}
};

static void 
mgga_c_vsxc_init(XC(func_type) *p)
{
  mgga_c_vsxc_params *params;

  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PW_MOD, XC_POLARIZED);

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_vsxc_params));
  params = (mgga_c_vsxc_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_C_VSXC:
    memcpy(params, &par_vsxc, sizeof(mgga_c_vsxc_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_vsxc\n");
    exit(1);
  }  
}


static void 
funcold(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static const FLOAT tmin = 0.5e-10;
  static const FLOAT sign[2] = {1.0, -1.0};

  mgga_c_vsxc_params *params;
  XC(lda_work_t) LDA[3];
  FLOAT opz, dd, g, dgdx, d2gdx2, h, dhdx, dhdt, aux, x_tot, dx_totdxs[2], ddddxs, ddddts;
  int is;

  assert(pt!=NULL && pt->params != NULL);
  params = (mgga_c_vsxc_params *)pt->params;

  /* first we get the parallel and perpendicular LDAS */
  XC(lda_stoll) (pt->func_aux[0], XC(lda_c_pw_func), r->dens, r->z, r->order, LDA);

  /* initialize to zero */
  r->f = 0.0;
  if(r->order >= 1){
    r->dfdrs = r->dfdz = r->dfdxs[0] = r->dfdxs[1] = r->dfdxt = 0.0;
    r->dfdus[0] = r->dfdus[1] = r->dfdts[0] = r->dfdts[1] = 0.0;
  }
  if(r->order >= 2){
    r->d2fdrs2 = r->d2fdrsz = r->d2fdrsxt = r->d2fdrsxs[0] = r->d2fdrsxs[1] = 0.0;
    r->d2fdz2 = r->d2fdzxt = r->d2fdzxs[0] = r->d2fdzxs[1] = r->d2fdxt2 = 0.0;
    r->d2fdxtxs[0] = r->d2fdxtxs[1] = r->d2fdxs2[0] = r->d2fdxs2[1] = r->d2fdxs2[2] = 0.0;
  }

  /* now we calculate the g functions for exchange and parallel correlation */
  for(is = 0; is < 2; is++){
    opz   = 1.0 + sign[is]*r->z;

    if(r->dens*opz < 2.0*pt->info->min_dens) continue;

    if(params->dss == NULL){
      h = dhdx = dhdt = 0.0;
    }else{
      XC(mgga_x_gvt4_func)(r->order, r->xs[is], 2.0*(r->ts[is] - K_FACTOR_C), 
			   params->alpha_ss, params->dss, &h, &dhdx, &dhdt);
    }

    g = dgdx = 0.0;

    dd = (r->ts[is] > tmin) ? 1.0 - r->xs[is]*r->xs[is]/(8.0*r->ts[is]) : 0.0;

    r->f += LDA[is].f*dd*(g + h);

    if(r->order < 1) continue;

    if(r->ts[is] > tmin){
      ddddxs = -2.0*r->xs[is]/(8.0*r->ts[is]);
      ddddts = r->xs[is]*r->xs[is]/(8.0*r->ts[is]*r->ts[is]);
    }else
      ddddxs = ddddts = 0.0;

    r->dfdrs     += LDA[is].dfdrs*dd*(g + h);
    r->dfdz      += LDA[is].dfdz *dd*(g + h);
    r->dfdxs[is] += LDA[is].f*(ddddxs*(g + h) + dd*(dgdx + dhdx));
    r->dfdts[is] += LDA[is].f*(ddddts*(g + h) + 2.0*dd*dhdt);
  }

  /* and now we add the opposite-spin contribution */
  aux   = r->xs[0]*r->xs[0] + r->xs[1]*r->xs[1];
  x_tot = SQRT(aux);

  if(params->dab == NULL){
    h = dhdx = dhdt = 0.0;
  }else{
    XC(mgga_x_gvt4_func)(r->order, x_tot, 2.0*(r->ts[0] + r->ts[1] - 2.0*K_FACTOR_C), 
			 params->alpha_ab, params->dab, &h, &dhdx, &dhdt);
  }

  g = dgdx = 0.0;

  r->f += LDA[2].f*(g + h);

  if(r->order < 1) return;

  dx_totdxs[0] = r->xs[0]/x_tot;
  dx_totdxs[1] = r->xs[1]/x_tot;

  r->dfdrs    += LDA[2].dfdrs*(g + h);
  r->dfdz     += LDA[2].dfdz *(g + h);
  r->dfdxs[0] += LDA[2].f*(dgdx + dhdx)*dx_totdxs[0];
  r->dfdxs[1] += LDA[2].f*(dgdx + dhdx)*dx_totdxs[1];
  r->dfdts[0] += LDA[2].f*dhdt*2.0;
  r->dfdts[1] += LDA[2].f*dhdt*2.0;
}

#include "maple2c/mgga_c_vsxc.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_vsxc) = {
  XC_MGGA_C_VSXC,
  XC_CORRELATION,
  "VSXC (correlation part)",
  XC_FAMILY_MGGA,
  {&xc_ref_VanVoorhis1998_400, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  0, NULL, NULL,
  mgga_c_vsxc_init,
  NULL,
  NULL, NULL,
  work_mgga_c,
};
