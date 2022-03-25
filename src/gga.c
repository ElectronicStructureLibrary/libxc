/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_gga.c"

/* Some useful formulas:

   sigma_st          = grad rho_s . grad rho_t
   zk                = energy density per unit particle

   vrho_s            = d zk / d rho_s
   vsigma_st         = d n*zk / d sigma_st

   v2rho2_st         = d^2 n*zk / d rho_s d rho_t
   v2rhosigma_svx    = d^2 n*zk / d rho_s d sigma_tv
   v2sigma2_stvx     = d^2 n*zk / d sigma_st d sigma_vx

   v3rho3_stv        = d^3 n*zk / d rho_s d rho_t d rho_v
   v3rho2sigma_stvx  = d^3 n*zk / d rho_s d rho_t d sigma_vx
   v3rhosigma2_svxyz = d^3 n*zk / d rho_s d sigma_vx d sigma_yz
   v3sigma3_stvxyz   = d^3 n*zk / d sigma_st d sigma_vx d sigma_yz

if nspin == 2
   rho(2)          = (u, d)
   sigma(3)        = (uu, ud, dd)

   vrho(2)         = (u, d)
   vsigma(3)       = (uu, ud, dd)

   v2rho2(3)       = (u_u, u_d, d_d)
   v2rhosigma(6)   = (u_uu, u_ud, u_dd, d_uu, d_ud, d_dd)
   v2sigma2(6)     = (uu_uu, uu_ud, uu_dd, ud_ud, ud_dd, dd_dd)

   v3rho3(4)       = (u_u_u, u_u_d, u_d_d, d_d_d)
   v3rho2sigma(9)  = (u_u_uu, u_u_ud, u_u_dd, u_d_uu, u_d_ud, u_d_dd, d_d_uu, d_d_ud, d_d_dd)
   v3rhosigma2(12) = (u_uu_uu, u_uu_ud, u_uu_dd, u_ud_ud, u_ud_dd, u_dd_dd, d_uu_uu, d_uu_ud, d_uu_dd, d_ud_ud, d_ud_dd, d_dd_dd)
   v3sigma(10)     = (uu_uu_uu, uu_uu_ud, uu_uu_dd, uu_ud_ud, uu_ud_dd, uu_dd_dd, ud_ud_ud, ud_ud_dd, ud_dd_dd, dd_dd_dd)

*/


void xc_evaluate_gga(const xc_func_type *func, int max_order,
                const xc_input_variables *in, xc_output_variables *out)
{
  int ii, check;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  /* turn off orders smaller than max_order */
  for(ii=max_order+1; ii <= XC_MAXIMUM_ORDER; ii++)
    orders[ii] = 0;

  /* check if all variables make sense */
  check = xc_input_variables_sanity_check(in, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    fprintf(stderr, "Field %s is not allocated\n", xc_input_variables_name[check]);
    exit(1);
  }
  
  check = xc_output_variables_sanity_check(out, orders, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    if(check >= 1000)
      fprintf(stderr, "Functional does not provide an implementation of the %d-th derivative\n", check-1000);
    else
      fprintf(stderr, "Field %s is not allocated\n", xc_output_variables_name[check]);
    exit(1);
  }
  
  xc_output_variables_initialize(out, in->np, func->nspin);
  
  /* call the GGA routines */
  if(func->info->gga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->gga->unpol[max_order] != NULL)
        func->info->gga->unpol[max_order](func, in->np, in->rho, in->sigma, out);
    }else{
      if(func->info->gga->pol[max_order] != NULL)
        func->info->gga->pol[max_order](func, in->np, in->rho, in->sigma, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, in, out);
}

/* old API */
void
xc_gga(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2,
       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3,
       double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4)
{
  int order = -1;

  if(zk     != NULL) order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;
  if(v3rho3 != NULL) order = 3;
  if(v4rho4 != NULL) order = 4;

  if(order < 0) return;

  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  out.v4rho4 = v4rho4; out.v4rho3sigma = v4rho3sigma; out.v4rho2sigma2 = v4rho2sigma2; out.v4rhosigma3 = v4rhosigma3; out.v4sigma4 = v4sigma4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, order, &in, &out);
}


/* specializations */
void
xc_gga_exc(const xc_func_type *p, size_t np, double *rho, double *sigma,
	    double *zk)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;
  
  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 0, &in, &out);
}

void
xc_gga_exc_vxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
		double *zk, double *vrho, double *vsigma)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;
  out.vrho = vrho; out.vsigma = vsigma;
  
  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 1, &in, &out);
}

void
xc_gga_exc_vxc_fxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
                    double *zk, double *vrho, double *vsigma,
                    double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  
  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 2, &in, &out);
}

void
xc_gga_vxc_fxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
                double *vrho, double *vsigma,
                double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  
  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 2, &in, &out);
}

void
xc_gga_exc_vxc_fxc_kxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
                        double *zk, double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                        double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;
  
  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 3, &in, &out);
}

void
xc_gga_vxc_fxc_kxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
                    double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                    double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho; out.vsigma = vsigma;
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 3, &in, &out);
}

void
xc_gga_vxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
           double *vrho, double *vsigma)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho; out.vsigma = vsigma;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 1, &in, &out);
}

void
xc_gga_fxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
           double *v2rho2, double *v2rhosigma, double *v2sigma2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v2rho2 = v2rho2; out.v2rhosigma = v2rhosigma; out.v2sigma2 = v2sigma2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 2, &in, &out);
}

void
xc_gga_kxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
           double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v3rho3 = v3rho3; out.v3rho2sigma = v3rho2sigma; out.v3rhosigma2 = v3rhosigma2; out.v3sigma3 = v3sigma3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 3, &in, &out);
}


void
xc_gga_lxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
           double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v4rho4 = v4rho4; out.v4rho3sigma = v4rho3sigma; out.v4rho2sigma2 = v4rho2sigma2; out.v4rhosigma3 = v4rhosigma3; out.v4sigma4 = v4sigma4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, sigma, NULL, NULL, NULL};

  xc_evaluate_gga(p, 4, &in, &out);
}
