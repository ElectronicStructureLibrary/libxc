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
set_ext_params_cpy_omega(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);

  /* This omega is only meant for internal use */
  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_NONE;
  p->hyb_coeff[0] = 0.0;
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams);
}

/*
   Copies parameters and sets the exact exchange coefficient, which
   should be the last parameter of the functional.
*/
void
set_ext_params_cpy_exx(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);

  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_FOCK;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[0] = 0.0;
}

/*
   Copies parameters and sets the HYB coefficients, which
   should be the three last parameters of the functional.
*/
void
set_ext_params_cpy_cam(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 3;
  copy_params(p, ext_params, nparams);

  assert(p->hyb_number_terms == 2);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams + 1);
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams + 2);

  p->hyb_type[1]  = XC_HYB_FOCK;
  p->hyb_coeff[1] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[1] = 0.0;
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
set_ext_params_cpy_cam_sr(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 2;
  copy_params(p, ext_params, nparams);

  assert(p->hyb_number_terms == 1);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = get_ext_param(p, ext_params, nparams);
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams + 1);
}

/* Long-range corrected functionals typically only have one parameter: the range separation parameter */
void
set_ext_params_cpy_lc(xc_func_type *p, const double *ext_params)
{
  int nparams;
  assert(p != NULL);
  nparams = p->info->ext_params.n - 1;
  copy_params(p, ext_params, nparams);

  assert(p->hyb_number_terms == 2);
  p->hyb_type[0]  = XC_HYB_ERF_SR;
  p->hyb_coeff[0] = -1.0;
  p->hyb_omega[0] = get_ext_param(p, ext_params, nparams);

  p->hyb_type[1]  = XC_HYB_FOCK;
  p->hyb_coeff[1] = 1.0;
  p->hyb_omega[1] = 0.0;
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

/* these are the dimension of the input and output
   arrays for spin unpolarized and polarized */
const xc_dimensions dimensions_unpolarized =
{
  /* inputs */
  1, 1, 1, 1, 1,    /* rho, sigma, lapl, tau, exx */
  /* order 0 */
  1,                /* zk */
  /* order 1 */
  1, 1, 1, 1, 1,    /* vrho, vsigma, vlapl, vtau, vexx */
  /* order 2 */
  1, 1, 1, 1, 1,    /* v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2rhoexx */
  1, 1, 1, 1,       /* v2sigma2, v2sigmalapl, v2sigmatau, v2sigmaexx */
  1, 1, 1,          /* v2lapl2, v2lapltau, v2laplexx */
  1, 1,             /* v2tau2, v2tauexx */
  1,                /* v2exx2 */
  /* order 3 */
  1, 1, 1, 1, 1,    /* v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rho2exx */
  1, 1, 1, 1,       /* v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rhosigmaexx */
  1, 1, 1,          /* v3rholapl2, v3rholapltau, v3rholaplexx */
  1, 1,             /* v3rhotau2, v3rhotauexx */
  1,                /* v3rhoexx2 */
  1, 1, 1, 1,       /* v3sigma3, v3sigma2lapl, v3sigma2tau, v3sigma2exx */
  1, 1, 1,          /* v3sigmalapl2, v3sigmalapltau, v3sigmalaplexx */
  1, 1,             /* v3sigmatau2, v3sigmatauexx */
  1,                /* v3sigmaexx2 */
  1, 1, 1,          /* v3lapl3, v3lapl2tau, v3lapl2exx */
  1, 1,             /* v3lapltau2, v3lapltauexx */
  1,                /* v3laplexx2 */
  1, 1, 1, 1,       /* v3tau3, v3tau2exx, v3tauexx2, v3exx3 */
  /* order 4 */
  1, 1, 1, 1, 1,    /* v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho3exx */
  1, 1, 1, 1,       /* v4rho2sigma2, v4rho2sigmalapl, v4rho2sigmatau, v4rho2sigmaexx */
  1, 1, 1,          /* v4rho2lapl2, v4rho2lapltau, v4rho2laplexx */
  1, 1,             /* v4rho2tau2, v4rho2tauexx */
  1,                /* v4rho2exx2 */
  1, 1, 1, 1,       /* v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau, v4rhosigma2exx */
  1, 1, 1,          /* v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmalaplexx */
  1, 1,             /* v4rhosigmatau2, v4rhosigmatauexx */
  1,                /* v4rhosigmaexx2 */
  1, 1, 1,          /* v4rhola1pl3, v4rholapl2tau, v4rholapl2exx */
  1, 1,             /* v4rholapltau2, v4rholapltauexx */
  1,                /* v4rholaplexx2 */
  1, 1, 1,          /* v4rhotau3, v4rhotau2exx, v4rhoexx3 */
  1, 1, 1, 1,       /* v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma3exx */
  1, 1, 1,          /* v4sigma2lapl2, v4sigma2lapltau, v4sigma2laplexx */
  1, 1,             /* v4sigma2tau2, v4sigma2tauexx */
  1,                /* v4sigma2exx2 */
  1, 1, 1,          /* v4sigmalapl3, v4sigmalapl2tau, v4sigmalapl2exx */
  1, 1,             /* v4sigmalapltau2, v4sigmalapltauexx */
  1,                /* v4sigmalaplexx2 */
  1, 1, 1, 1,       /* v4sigmatau3, v4sigmatau2exx, v4sigmatauexx2, v4sigmaexx3 */
  1, 1, 1,          /* v4lapl4, v4lapl3tau, v4lapl3exx */
  1, 1, 1,          /* v4lapl2tau2, v4lapl2tauexx, v4lapl2exx2 */
  1, 1, 1, 1,       /* v4lapltau3, v4lapltau2exx, v4lapltauexx2, v4laplexx3 */
  1, 1, 1, 1        /* v4tau4, v4tau3exx, v4tauexx3, v4exx4 */
};

const xc_dimensions dimensions_polarized =
{
  /* inputs */
  2, 3, 2, 2, 2,    /* rho, sigma, lapl, tau, exx */
  /* order 0 */
  1,                /* zk */
  /* order 1 */
  2, 3, 2, 2, 2,    /* vrho, vsigma, vlapl, vtau, vexx */
  /* order 2 */
  3, 6, 4, 4, 4,    /* v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2rhoexx */
  6, 6, 6, 6,       /* v2sigma2, v2sigmalapl, v2sigmatau, v2sigmaexx */
  3, 4, 4,          /* v2lapl2, v2lapltau, v2laplexx */
  3, 4,             /* v2tau2, v2tauexx */
  3,                /* v2exx2 */
  /* order 3 */
  4, 9, 6, 6, 6,    /* v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rho2exx */
  12, 12, 12, 12,   /* v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rhosigmaexx */
  6, 8, 8,          /* v3rholapl2, v3rholapltau, v3rholaplexx */
  6, 8,             /* v3rhotau2, v3rhotauexx */
  6,                /* v3rhoexx2 */
  10, 12, 12, 12,   /* v3sigma3, v3sigma2lapl, v3sigma2tau, v3sigma2exx */
  9, 12, 12,        /* v3sigmalapl2, v3sigmalapltau, v3sigmalaplexx */
  9, 12,            /* v3sigmatau2, v3sigmatauexx */
  9,                /* v3sigmaexx2 */
  4, 6, 6,          /* v3lapl3, v3lapl2tau, v3lapl2exx */
  6, 8,             /* v3lapltau2, v3lapltauexx */
  6,                /* v3laplexx2 */
  4, 6, 6, 4,       /* v3tau3, v3tau2exx, v3tauexx2, v3exx3 */
  /* order 4 */
  5, 12, 8, 8, 8,   /* v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho3exx */
  18, 18, 18, 18,   /* v4rho2sigma2, v4rho2sigmalapl, v4rho2sigmatau, v4rho2sigmaexx */
  9, 12, 12,        /* v4rho2lapl2, v4rho2lapltau, v4rho2laplexx */
  9, 12,            /* v4rho2tau2, v4rho2tauexx */
  9,                /* v4rho2exx2 */
  20, 36, 36, 36,   /* v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau, v4rhosigma2exx */
  18, 24, 24,       /* v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmalaplexx */
  18, 24,           /* v4rhosigmatau2, v4rhosigmatauexx */
  18,               /* v4rhosigmaexx2 */
  8, 12, 12,        /* v4rholapl3, v4rholapl2tau, v4rholapl2exx */
  12, 16,           /* v4rholapltau2, v4rholapltauexx */
  12,               /* v4rholaplexx2 */
  8, 12, 8,         /* v4rhotau3, v4rhotau2exx, v4rhoexx3 */
  15, 20, 20, 20,   /* v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma3exx */
  18, 24, 24,       /* v4sigma2lapl2, v4sigma2lapltau, v4sigma2laplexx */
  18, 24,           /* v4sigma2tau2, v4sigma2tauexx */
  18,               /* v4sigma2exx2 */
  12, 18, 18,       /* v4sigmalapl3, v4sigmalapl2tau, v4sigmalapl2exx */
  18, 24,           /* v4sigmalapltau2, v4sigmalapltauexx */
  18,               /* v4sigmalaplexx2 */
  12, 18, 18, 12,   /* v4sigmatau3, v4sigmatau2exx, v4sigmatauexx2, v4sigmaexx3 */
  5, 8, 8,          /* v4lapl4, v4lapl3tau, v4lapl3exx */
  9, 12, 9,         /* v4lapl2tau2, v4lapl2tauexx, v4lapl2exx2 */
  8, 12, 12, 8,     /* v4lapltau3, v4lapltau2exx, v4lapltauexx2, v4laplexx3 */
  5, 8, 8, 5        /* v4tau4, v4tau3exx, v4tauexx3, v4exx4 */
};

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
