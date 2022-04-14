/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_lda.c"

/* old API */
void
xc_lda(const xc_func_type *p, size_t np,
       double *rho,
       double *zk, double *vrho, double *v2rho2, double *v3rho3, double *v4rho4)
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
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;
  out.v4rho4 = v4rho4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};
  
  xc_evaluate_func(p, order, &in, &out);
}


/* specializations */
void
xc_lda_exc(const xc_func_type *p, size_t np, double *rho, double *zk)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 0, &in, &out);
}

void
xc_lda_exc_vxc(const xc_func_type *p, size_t np, double *rho, double *zk, double *vrho)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;
  out.vrho = vrho;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 1, &in, &out);
}

void
xc_lda_exc_vxc_fxc(const xc_func_type *p, size_t np, double *rho, double *zk, double *vrho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 2, &in, &out);
}

void
xc_lda_vxc_fxc(const xc_func_type *p, size_t np, double *rho, double *vrho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 2, &in, &out);
}

void
xc_lda_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 3, &in, &out);
}

void
xc_lda_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 3, &in, &out);
}

void
xc_lda_vxc(const xc_func_type *p, size_t np, double *rho, double *vrho)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 1, &in, &out);
}

void
xc_lda_fxc(const xc_func_type *p, size_t np, double *rho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v2rho2 = v2rho2;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};
  
  xc_evaluate_func(p, 2, &in, &out);
}

void
xc_lda_kxc(const xc_func_type *p, size_t np, double *rho, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v3rho3 = v3rho3;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 3, &in, &out);
}

void
xc_lda_lxc(const xc_func_type *p, size_t np, double *rho, double *v4rho4)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v4rho4 = v4rho4;

  const xc_input_variables_dimensions *idim = input_variables_dimensions_get(p->nspin);
  const xc_input_variables in =
    {np, idim, rho, NULL, NULL, NULL, NULL};

  xc_evaluate_func(p, 4, &in, &out);
}
