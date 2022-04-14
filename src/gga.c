/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_gga.c"

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

  xc_evaluate_func(p, order, &in, &out);
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

  xc_evaluate_func(p, 0, &in, &out);
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

  xc_evaluate_func(p, 1, &in, &out);
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

  xc_evaluate_func(p, 2, &in, &out);
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

  xc_evaluate_func(p, 2, &in, &out);
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

  xc_evaluate_func(p, 3, &in, &out);
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

  xc_evaluate_func(p, 3, &in, &out);
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

  xc_evaluate_func(p, 1, &in, &out);
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

  xc_evaluate_func(p, 2, &in, &out);
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

  xc_evaluate_func(p, 3, &in, &out);
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

  xc_evaluate_func(p, 4, &in, &out);
}
