/*
 Copyright (C) 2006-2009 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_MGGA_X_BJ06         207 /* Becke & Johnson correction to Becke-Roussel 89  */
#define XC_MGGA_X_TB09         208 /* Tran & Blaha correction to Becke & Johnson  */
#define XC_MGGA_X_RPP09        209 /* Rasanen, Pittalis, and Proetto correction to Becke & Johnson  */

typedef struct xc_mgga_work_x_t {
  int   order; /* to which order should I return the derivatives */
  double rs, zeta, x, t, u;

  double f;                                   /* enhancement factor       */
  double dfdrs, dfdx, dfdt, dfdu;             /* first derivatives of f  */
  double d2fdrs2, d2fdx2, d2fdt2, d2fdu2;     /* second derivatives of zk */
  double d2fdrsx, d2fdrst, d2fdrsu, d2fdxt, d2fdxu, d2fdtu;
} xc_mgga_work_x_t;

typedef struct{
  double c;
} mgga_x_tb09_params;

static double br89_gamma = 0.8;
static double b00_at     = 0.928;

static void
mgga_x_tb09_init(xc_func_type *p)
{
  mgga_x_tb09_params *params;

  p->params = malloc(sizeof(mgga_x_tb09_params));
  params = (mgga_x_tb09_params *)p->params;

  params->c = 0;

  switch(p->info->number){
  case XC_MGGA_X_BJ06:
    params->c = 1.0;
    break;
  case XC_MGGA_X_TB09:
    /* the value of c should be passed by the calling code */
    break;
  case XC_MGGA_X_RPP09:
    params->c = 1.0;
    break;
  }
}

static void
func(const xc_func_type *pt, xc_mgga_work_x_t *r)
{
  double Q, br_x, v_BR, dv_BRdbx, d2v_BRdbx2, dxdQ, d2xdQ2, ff, dffdx, d2ffdx2;
  double cnst, c_TB09, c_HEG, exp1, exp2, gamma, fw, dfwdt, min_Q;

  min_Q = 5.0e-13;

  gamma = (pt->info->number == XC_MGGA_X_B00 || pt->info->number == XC_MGGA_X_BR89_1) ? 1.0 : br89_gamma;

  Q = (r->u - 4.0*gamma*r->t + 0.5*gamma*r->x*r->x)/6.0;
  if(fabs(Q) < min_Q) Q = (Q < 0) ? -min_Q : min_Q;

  br_x = xc_mgga_x_br89_get_x(Q);

  cnst = -2.0*CBRT(M_PI)/X_FACTOR_C;
  exp1 = exp(br_x/3.0);
  exp2 = exp(-br_x);

  v_BR = (fabs(br_x) > pt->dens_threshold) ?
    exp1*(1.0 - exp2*(1.0 + br_x/2.0))/br_x :
    1.0/2.0 + br_x/6.0 - br_x*br_x/18.0;

  v_BR *= cnst;

  if(pt->info->number == XC_MGGA_X_BR89 || pt->info->number == XC_MGGA_X_BR89_1){
    /* we have also to include the factor 1/2 from Eq. (9) */
    r->f = - v_BR / 2.0;

  }else{ /* XC_MGGA_X_BJ06 & XC_MGGA_X_TB09 */
    r->f = 0.0;
  }

  if(r->order < 1) return;

  if(pt->info->number == XC_MGGA_X_BR89 || pt->info->number == XC_MGGA_X_BR89_1 || r->order > 1){
    dv_BRdbx = (fabs(br_x) > pt->dens_threshold) ?
      (3.0 + br_x*(br_x + 2.0) + (br_x - 3.0)/exp2) / (3.0*exp1*exp1*br_x*br_x) :
      1.0/6.0 - br_x/9.0;
    dv_BRdbx *= cnst;

    ff    = br_x*exp(-2.0/3.0*br_x)/(br_x - 2);
    dffdx = ff*(-2.0/3.0 + 1.0/br_x - 1.0/(br_x - 2.0));
    dxdQ  = -ff/(Q*dffdx);
  }

  if(pt->info->number == XC_MGGA_X_BR89 || pt->info->number == XC_MGGA_X_BR89_1){
    r->dfdx = -r->x*gamma*dv_BRdbx*dxdQ/12.0;
    r->dfdt =   4.0*gamma*dv_BRdbx*dxdQ/12.0;
    r->dfdu =            -dv_BRdbx*dxdQ/12.0;

  }else{
    assert(pt->params != NULL);
    c_TB09 = ((mgga_x_tb09_params *) (pt->params))->c;

    r->dfdrs = -c_TB09*v_BR;

    c_HEG  = (3.0*c_TB09 - 2.0)*sqrt(5.0/12.0)/(X_FACTOR_C*M_PI);

    if(pt->info->number == XC_MGGA_X_BJ06 || pt->info->number == XC_MGGA_X_TB09)
      r->dfdrs -= c_HEG*sqrt(2.0*r->t);
    else /* XC_MGGA_X_RPP09 */
      r->dfdrs -= c_HEG*sqrt(max(2.0*r->t - r->x*r->x/4.0, 0.0));

    r->dfdrs /= -r->rs; /* due to the definition of dfdrs */
  }

  if(r->order < 2) return;

  if(pt->info->number == XC_MGGA_X_BR89 || pt->info->number == XC_MGGA_X_BR89_1 || r->order > 2){
    d2v_BRdbx2 = (fabs(br_x) > pt->dens_threshold) ?
      ((18.0 + (br_x - 6.0)*br_x)/exp2 - 2.0*(9.0 + br_x*(6.0 + br_x*(br_x + 2.0))))
      / (9.0*exp1*exp1*br_x*br_x*br_x) :
      -1.0/9.0;
    d2v_BRdbx2 *= cnst;

    d2ffdx2 = dffdx*dffdx/ff + ff*(-1.0/(br_x*br_x) + 1.0/((br_x - 2.0)*(br_x - 2.0)));
    d2xdQ2 = -(2.0*dxdQ/Q + d2ffdx2*dxdQ*dxdQ/dffdx);
  }

  if(pt->info->number == XC_MGGA_X_BR89 || pt->info->number == XC_MGGA_X_BR89_1){
    double aux1 = d2v_BRdbx2*dxdQ*dxdQ + dv_BRdbx*d2xdQ2;

    r->d2fdx2 = -(aux1*gamma*r->x*r->x/6.0 + dv_BRdbx*dxdQ)*gamma/12.0;
    r->d2fdxt =  aux1*gamma*gamma*r->x/18.0;
    r->d2fdxu = -aux1*gamma*r->x/72.0;
    r->d2fdt2 = -aux1*2.0*gamma*gamma/9.0;
    r->d2fdtu =  aux1*gamma/18.0;
    r->d2fdu2 = -aux1/72.0;
  }else{

  }

}

#include "work_mgga_x.c"

const xc_func_info_type xc_func_info_mgga_x_bj06 = {
  XC_MGGA_X_BJ06,
  XC_EXCHANGE,
  "Becke & Johnson 06",
  XC_FAMILY_MGGA,
  {&xc_ref_Becke2006_221101, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_VXC,
  1e-23,
  0, NULL, NULL,
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga_x,
};

static const func_params_type ext_params[] = {
  {"c", 1.0, "This parameter involves an average over the unit cell and must be calculated by the calling program."},
};

static void
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_tb09_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_tb09_params *) (p->params);

  params->c = get_ext_param(p->info->ext_params, ext_params, 0);
}

const xc_func_info_type xc_func_info_mgga_x_tb09 = {
  XC_MGGA_X_TB09,
  XC_EXCHANGE,
  "Tran & Blaha 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Tran2009_226401, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_VXC,
  1.0e-23,
  1, ext_params, set_ext_params,
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga_x,
};

const xc_func_info_type xc_func_info_mgga_x_rpp09 = {
  XC_MGGA_X_RPP09,
  XC_EXCHANGE,
  "Rasanen, Pittalis & Proetto 09",
  XC_FAMILY_MGGA,
  {&xc_ref_Rasanen2010_044112, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_I_HAVE_VXC,
  1e-23,
  0, NULL, NULL,
  mgga_x_tb09_init, NULL,
  NULL, NULL, work_mgga_x,
};


