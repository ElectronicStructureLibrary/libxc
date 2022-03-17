/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_lda.c"

void
xc_lda_sanity_check(const xc_func_info_type *info, int order, xc_output_variables *out)
{
  /* sanity check */
  if(order < 0 || order > 4){
    fprintf(stderr, "Order of derivatives '%d' not implemented\n",
	    order);
    exit(1);
  }

  if(out->zk != NULL && !(info->flags & XC_FLAGS_HAVE_EXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of Exc\n",
	    info->name);
    exit(1);
  }

  if(out->vrho != NULL && !(info->flags & XC_FLAGS_HAVE_VXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of vxc\n",
	    info->name);
    exit(1);
  }

  if(out->v2rho2 != NULL && !(info->flags & XC_FLAGS_HAVE_FXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of fxc\n",
	    info->name);
    exit(1);
  }

  if(out->v3rho3 != NULL && !(info->flags & XC_FLAGS_HAVE_KXC)){
    fprintf(stderr, "Functional '%s' does not provide an implementation of kxc\n",
	    info->name);
    exit(1);
  }
}


/* get the lda functional */
void
xc_evaluate_lda(const xc_func_type *func, int order, size_t np, const double *rho,
       xc_output_variables *out)
{
  xc_lda_sanity_check(func->info, order, out);
  xc_output_variables_initialize(out, np, func->nspin);

  /* call the LDA routines */
  if(func->info->lda != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->lda->unpol[order] != NULL)
        func->info->lda->unpol[order](func, np, rho, out);
    }else{
      if(func->info->lda->pol[order] != NULL)
        func->info->lda->pol[order](func, np, rho, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, np, rho, NULL, NULL, NULL, NULL, out);
}

/* old API */
void
xc_lda(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3, double *v4rho4)
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

  xc_evaluate_lda(p, order, np, rho, &out);
}


/* specializations */
void
xc_lda_exc(const xc_func_type *p, size_t np, const double *rho, double *zk)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;

  xc_evaluate_lda(p, 0, np, rho, &out);
}

void
xc_lda_exc_vxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk   = zk;
  out.vrho = vrho;

  xc_evaluate_lda(p, 1, np, rho, &out);
}

void
xc_lda_exc_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  xc_evaluate_lda(p, 2, np, rho, &out);
}

void
xc_lda_vxc_fxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;

  xc_evaluate_lda(p, 2, np, rho, &out);
}

void
xc_lda_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.zk     = zk;
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  xc_evaluate_lda(p, 3, np, rho, &out);
}

void
xc_lda_vxc_fxc_kxc(const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;
  out.v2rho2 = v2rho2;
  out.v3rho3 = v3rho3;

  xc_evaluate_lda(p, 3, np, rho, &out);
}

void
xc_lda_vxc(const xc_func_type *p, size_t np, const double *rho, double *vrho)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.vrho   = vrho;

  xc_evaluate_lda(p, 1, np, rho, &out);
}

void
xc_lda_fxc(const xc_func_type *p, size_t np, const double *rho, double *v2rho2)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v2rho2 = v2rho2;

  xc_evaluate_lda(p, 2, np, rho, &out);
}

void
xc_lda_kxc(const xc_func_type *p, size_t np, const double *rho, double *v3rho3)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v3rho3 = v3rho3;

  xc_evaluate_lda(p, 3, np, rho, &out);
}

void
xc_lda_lxc(const xc_func_type *p, size_t np, const double *rho, double *v4rho4)
{
  xc_output_variables out;
  libxc_memset(&out, 0, sizeof(xc_output_variables));
  out.v4rho4 = v4rho4;

  xc_evaluate_lda(p, 4, np, rho, &out);
}
