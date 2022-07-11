/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2015-2021 Susi Lehtola
               2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"


/* this function converts the spin-density into total density and
	 relative magnetization */
/* inline */ GPU_FUNCTION void
xc_rho2dzeta(int nspin, const double *rho, double *d, double *zeta)
{
  if(nspin==XC_UNPOLARIZED){
    *d    = m_max(rho[0], 0.0);
    *zeta = 0.0;
  }else{
    *d = rho[0] + rho[1];
    if(*d > 0.0){
      *zeta = (rho[0] - rho[1])/(*d);
      *zeta = m_min(*zeta,  1.0);
      *zeta = m_max(*zeta, -1.0);
    }else{
      *d    = 0.0;
      *zeta = 0.0;
    }
  }
}

const char *get_kind(const xc_func_type *func) {
  switch(func->info->kind) {
   case(XC_EXCHANGE):
      return "XC_EXCHANGE";

    case(XC_CORRELATION):
      return "XC_CORRELATION";

    case(XC_EXCHANGE_CORRELATION):
      return "XC_EXCHANGE_CORRELATION";

    case(XC_KINETIC):
      return "XC_KINETIC";

    default:
      printf("Internal error in get_kind.\n");
      return "";
  }
}

const char *get_family(const xc_func_type *func) {
  switch(func->info->family) {
    case(XC_FAMILY_UNKNOWN):
      return "XC_FAMILY_UNKNOWN";

    case(XC_FAMILY_LDA):
      return "XC_FAMILY_LDA";

    case(XC_FAMILY_GGA):
      return "XC_FAMILY_GGA";

    case(XC_FAMILY_MGGA):
      return "XC_FAMILY_MGGA";

    case(XC_FAMILY_LCA):
      return "XC_FAMILY_LCA";

    case(XC_FAMILY_OEP):
      return "XC_FAMILY_OEP";

    default:
      printf("Internal error in get_family.\n");
      return "";
  }
}

/* this function checks if it should use the default or
   the user assigned value for an external parameter */
double
get_ext_param(const xc_func_type *func, const double *values, int index)
{
  assert(index >= 0 && index < func->info->ext_params.n);
  return func->ext_params[index];
}

/* Copy n parameters, assumes that p->params is just a series of doubles
   so it can be accessed as a array, and and copies
   ext_params to this. */
static void copy_params(xc_func_type *p, const double *ext_params, int nparams) {
  double *params;
  int ii;

  assert(nparams >= 0);
  if(nparams) {
    /* Some functionals only set the hybrid parameters which require no extra storage */
    assert(p->params != NULL);
    params = (double *) (p->params);
    for(ii=0; ii<nparams; ii++)
      params[ii] = get_ext_param(p, ext_params, ii);
  }
}

/* Just copy the parameters */
void
set_ext_params_cpy(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n;
  copy_params(p, ext_params, nparams);
}

/*
   Copies parameters and sets the screening parameter, which should be
   the last parameter of the functional.
*/
void
set_ext_params_omega(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;

  /* This omega is only meant for internal use */
  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_NONE;
  p->hyb_coeff[0] = 0.0;
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams);
}

void
set_ext_params_cpy_omega(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);
  set_ext_params_omega(p, ext_params);
}

/*
   Copies parameters and sets the exact exchange coefficient, which
   should be the last parameter of the functional.
*/
void
set_ext_params_exx(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;

  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_FOCK;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[0] = 0.0;
}

void
set_ext_params_cpy_exx(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);
  set_ext_params_exx(p, ext_params);
}

/*
   Copies parameters and sets the HYB coefficients, which
   should be the three last parameters of the functional.
*/
void
set_ext_params_cam(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 3;

  assert(p->hyb_number_terms == 2);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams + 1);
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams + 2);

  p->hyb_type[1]  = XC_HYB_FOCK;
  p->hyb_coeff[1] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[1] = 0.0;
}

void
set_ext_params_cpy_cam(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 3;
  copy_params(p, ext_params, nparams);
  set_ext_params_cam(p, ext_params);
}

void
set_ext_params_camy(xc_func_type *p, const double *ext_params)
{
  set_ext_params_cam(p, ext_params);
  p->hyb_type[0]  = XC_HYB_YUKAWA_SR;
}

void
set_ext_params_cpy_camy(xc_func_type *p, const double *ext_params)
{
  set_ext_params_cpy_cam(p, ext_params);
  p->hyb_type[0]  = XC_HYB_YUKAWA_SR;
}

/*
  Short-range-only version
*/
void
set_ext_params_cam_sr(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 2;

  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams + 1);
}

void
set_ext_params_cpy_cam_sr(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 2;
  copy_params(p, ext_params, nparams);
  set_ext_params_cam_sr(p, ext_params);
}

/* Long-range corrected functionals typically only have one parameter: the range separation parameter */
void
set_ext_params_lc(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;

  assert(p->hyb_number_terms == 2);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = -1.0;
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams);

  p->hyb_type[1]  = XC_HYB_FOCK;
  p->hyb_coeff[1] = 1.0;
  p->hyb_omega[1] = 0.0;
}

void
set_ext_params_cpy_lc(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);
  set_ext_params_lc(p, ext_params);
}

void
set_ext_params_lcy(xc_func_type *p, const double *ext_params)
{
  set_ext_params_lc(p, ext_params);
  p->hyb_type[0]  = XC_HYB_YUKAWA_SR;
}

void
set_ext_params_cpy_lcy(xc_func_type *p, const double *ext_params)
{
  set_ext_params_cpy_lc(p, ext_params);
  p->hyb_type[0]  = XC_HYB_YUKAWA_SR;
}

/* Free pointer */
void
libxc_free(void *ptr)
{
#ifdef HAVE_CUDA
  cudaFree(ptr);
#else
  free(ptr);
#endif
}

/** Computes nderiv derivatives of B-spline Nip(u)

    The algorithm follows the textbook presentation in the NURBS
    book, 2nd edition, by Les Piegl and Wayne Tiller.

    Input variables:
    - i: function index
    - p: spline order
    - u: argument
    - nderiv: number of derivatives to calculate (zero for just the function itself)
    - U: knots
    Output variables:
    - ders: array [Nip(u), Nip'(u), Nip''(u), ..., Nip^(nderiv)(u)]
*/
GPU_FUNCTION void
xc_bspline(int i, int p, double u, int nderiv, const double *U, double *ders) {

  /* Initialize output array */
  libxc_memset(ders, 0, (nderiv+1)*sizeof(double));

  /* Check locality of support */
  if(u < U[i] || u >= U[i+p+1]) {
    return;
  }

  /* Arrays need static sizes for stack allocation */
#define PMAX 8
  assert(p<PMAX);

  /* Array of normalized B splines, use dense storage for simpler code */
  double N[PMAX][PMAX];
  libxc_memset(N, 0, PMAX*PMAX*sizeof(double));

  /* Initialize zeroth-degree functions: piecewise constants */
  for(int j=0; j<=p; j++) {
    N[0][j] = (u >= U[i+j] && u < U[i+j+1]) ? 1.0 : 0.0;
  }

  /* Fill out table of B splines */
  for(int k=1; k<=p; k++) {
    double saved = (N[k-1][0] == 0.0) ? 0.0 : ((u-U[i])*N[k-1][0])/(U[i+k]-U[i]);

    for(int j=0; j<=p-k; j++) {
      double Ul = U[i+j+1];
      double Ur = U[i+j+k+1];
      if(N[k-1][j+1] == 0.0) {
        N[k][j] = saved;
        saved = 0.0;
      } else {
        double temp = N[k-1][j+1] / (Ur-Ul);
        N[k][j] = saved + (Ur-u)*temp;
        saved = (u-Ul)*temp;
      }
    }
  }

  /* Function value */
  ders[0] = N[p][0];
  if(nderiv==0)
    return;

  /* Helper memory */
  assert(nderiv<=4);
  double ND[5]; /* dimension nderiv+1 */
  int maxk = (nderiv < p) ? nderiv : p;

  /* Compute derivatives */
  for(int k=1; k<=maxk; k++) {
    /* Load appropriate column */
    libxc_memset(ND, 0, (nderiv+1)*sizeof(double));
    for(int j=0; j<=k; j++)
      ND[j] = N[p-k][j];

    /* Compute table */
    for(int jj=1; jj<=k; jj++) {
      double saved = (ND[0] == 0.0) ? 0.0 : ND[0]/(U[i+p-k+jj]-U[i]);

      for(int j=0; j<=k-jj; j++) {
        double Ul = U[i+j+1];
        /* the -k term is missing in the book */
        double Ur = U[i+j+p-k+jj+1];
        if(ND[j+1] == 0.0) {
          ND[j] = (p-k+jj)*saved;
          saved = 0.0;
        } else {
          double temp = ND[j+1]/(Ur-Ul);
          ND[j] = (p-k+jj)*(saved-temp);
          saved = temp;
        }
      }
    }
    /* k:th derivative is */
    ders[k] = ND[0];
  }
}
