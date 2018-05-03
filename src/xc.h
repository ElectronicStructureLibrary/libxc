/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef _XC_H
#define _XC_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xc_version.h>

#define XC_UNPOLARIZED          1
#define XC_POLARIZED            2

#define XC_NON_RELATIVISTIC     0
#define XC_RELATIVISTIC         1

#define XC_EXCHANGE             0
#define XC_CORRELATION          1
#define XC_EXCHANGE_CORRELATION 2
#define XC_KINETIC              3

#define XC_FAMILY_UNKNOWN      -1
#define XC_FAMILY_LDA           1
#define XC_FAMILY_GGA           2
#define XC_FAMILY_MGGA          4
#define XC_FAMILY_LCA           8
#define XC_FAMILY_OEP          16
#define XC_FAMILY_HYB_GGA      32
#define XC_FAMILY_HYB_MGGA     64

/* flags that can be used in info.flags. Don't reorder these since it
   will break the ABI of the library. */
#define XC_FLAGS_HAVE_EXC         (1 <<  0) /*     1 */
#define XC_FLAGS_HAVE_VXC         (1 <<  1) /*     2 */
#define XC_FLAGS_HAVE_FXC         (1 <<  2) /*     4 */
#define XC_FLAGS_HAVE_KXC         (1 <<  3) /*     8 */
#define XC_FLAGS_HAVE_LXC         (1 <<  4) /*    16 */
#define XC_FLAGS_1D               (1 <<  5) /*    32 */
#define XC_FLAGS_2D               (1 <<  6) /*    64 */
#define XC_FLAGS_3D               (1 <<  7) /*   128 */
#define XC_FLAGS_HYB_CAM          (1 <<  8) /*   256 */
#define XC_FLAGS_HYB_CAMY         (1 <<  9) /*   512 */
#define XC_FLAGS_VV10             (1 << 10) /*  1024 */
#define XC_FLAGS_HYB_LC           (1 << 11) /*  2048 */
#define XC_FLAGS_HYB_LCY          (1 << 12) /*  4096 */
#define XC_FLAGS_STABLE           (1 << 13) /*  8192 */
#define XC_FLAGS_DEVELOPMENT      (1 << 14) /* 16384 */
#define XC_FLAGS_NEEDS_LAPLACIAN  (1 << 15) /* 32768 */

#define XC_TAU_EXPLICIT         0
#define XC_TAU_EXPANSION        1

#define XC_MAX_REFERENCES       5

void xc_version(int *major, int *minor, int *micro);
const char *xc_version_string();

struct xc_func_type;

typedef struct{
  char *ref, *doi, *bibtex;
} func_reference_type;

char const *xc_func_reference_get_ref(const func_reference_type *reference);
char const *xc_func_reference_get_doi(const func_reference_type *reference);
char const *xc_func_reference_get_bibtex(const func_reference_type *reference);

typedef struct{
  double value;
  char *description;
} func_params_type;

typedef struct{
  int   number;   /* identifier number */
  int   kind;     /* XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC */

  char *name;     /* name of the functional, e.g. "PBE" */
  int   family;   /* type of the functional, e.g. XC_FAMILY_GGA */
  func_reference_type *refs[XC_MAX_REFERENCES];  /* index of the references */

  int   flags;    /* see above for a list of possible flags */

  double dens_threshold;

  /* this allows to have external parameters in the functional */
  int n_ext_params;
  const func_params_type *ext_params;
  void (*set_ext_params)(struct xc_func_type *p, const double *ext_params);

  void (*init)(struct xc_func_type *p);
  void (*end) (struct xc_func_type *p);
  void (*lda) (const struct xc_func_type *p, int np,
	       const double *rho,
	       double *zk, double *vrho, double *v2rho2, double *v3rho3);
  void (*gga) (const struct xc_func_type *p, int np,
	       const double *rho, const double *sigma,
	       double *zk, double *vrho, double *vsigma,
	       double *v2rho2, double *v2rhosigma, double *v2sigma2,
	       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
  void (*mgga)(const struct xc_func_type *p, int np,
	       const double *rho, const double *sigma, const double *lapl_rho, const double *tau,
	       double *zk, double *vrho, double *vsigma, double *vlapl_rho, double *vtau,
	       double *v2rho2, double *v2sigma2, double *v2lapl2, double *v2tau2,
	       double *v2rhosigma, double *v2rholapl, double *v2rhotau,
	       double *v2sigmalapl, double *v2sigmatau, double *v2lapltau);
} xc_func_info_type;

/* for API compability with older versions of libxc */
#define XC(func) xc_ ## func

int xc_func_info_get_number(const xc_func_info_type *info);
int xc_func_info_get_kind(const xc_func_info_type *info);
char const *xc_func_info_get_name(const xc_func_info_type *info);
int xc_func_info_get_family(const xc_func_info_type *info);
int xc_func_info_get_flags(const xc_func_info_type *info);
const func_reference_type *xc_func_info_get_references(const xc_func_info_type *info, int number);
int xc_func_info_get_n_ext_params(xc_func_info_type *info);
char const *xc_func_info_get_ext_params_description(xc_func_info_type *info, int number);
double xc_func_info_get_ext_params_default_value(xc_func_info_type *info, int number);

struct xc_func_type{
  const xc_func_info_type *info;       /* all the information concerning this functional */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */

  int n_func_aux;                       /* how many auxiliary functions we need */
  struct xc_func_type **func_aux;      /* most GGAs are based on a LDA or other GGAs  */
  double *mix_coef;                      /* coefficients for the mixing */

  /**
     Parameters for range-separated hybrids
     cam_omega: the range separation constant
     cam_alpha: fraction of full Hartree-Fock exchange, used both for
                usual hybrids as well as range-separated ones
     cam_beta:  fraction of short-range only(!) exchange in
                range-separated hybrids

     N.B. Different conventions for alpha and beta can be found in
     literature. In the convention used in libxc, at short range the
     fraction of exact exchange is cam_alpha+cam_beta, while at long
     range it is cam_alpha.
  */
  double cam_omega, cam_alpha, cam_beta;

  double nlc_b;                          /* Non-local correlation, b parameter */
  double nlc_C;                          /* Non-local correlation, C parameter */

  int n_rho, n_sigma, n_tau, n_lapl;    /* spin dimensions of the arrays */
  int n_zk;

  int n_vrho, n_vsigma, n_vtau, n_vlapl;

  int n_v2rho2, n_v2sigma2, n_v2tau2, n_v2lapl2,
    n_v2rhosigma, n_v2rhotau, n_v2rholapl,
    n_v2sigmatau, n_v2sigmalapl, n_v2lapltau;

  int n_v3rho3, n_v3rho2sigma, n_v3rhosigma2, n_v3sigma3;

  void *params;                         /* this allows us to fix parameters in the functional */
  double dens_threshold;
};

typedef struct xc_func_type xc_func_type;

/* functionals */
int   xc_functional_get_number(const char *name);
char *xc_functional_get_name(int number);
int   xc_family_from_id(int id, int *family, int *number);
int   xc_number_of_functionals();
int   xc_maximum_name_length();
void xc_available_functional_numbers(int *list);
void xc_available_functional_names(char **list);

  xc_func_type *xc_func_alloc();
int   xc_func_init(xc_func_type *p, int functional, int nspin);
void  xc_func_end(xc_func_type *p);
void  xc_func_free(xc_func_type *p);
const xc_func_info_type *xc_func_get_info(const xc_func_type *p);
void xc_func_set_dens_threshold(xc_func_type *p, double dens_threshold);
void  xc_func_set_ext_params(xc_func_type *p, double *ext_params);

#include "xc_funcs.h"
#include "xc_funcs_removed.h"

void xc_lda        (const xc_func_type *p, int np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3);
void xc_lda_exc    (const xc_func_type *p, int np, const double *rho, double *zk);
void xc_lda_exc_vxc(const xc_func_type *p, int np, const double *rho, double *zk, double *vrho);
void xc_lda_vxc    (const xc_func_type *p, int np, const double *rho, double *vrho);
void xc_lda_fxc    (const xc_func_type *p, int np, const double *rho, double *v2rho2);
void xc_lda_kxc    (const xc_func_type *p, int np, const double *rho, double *v3rho3);

void xc_gga     (const xc_func_type *p, int np, const double *rho, const double *sigma,
		  double *zk, double *vrho, double *vsigma,
		  double *v2rho2, double *v2rhosigma, double *v2sigma2,
		  double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_gga_exc(const xc_func_type *p, int np, const double *rho, const double *sigma,
		 double *zk);
void xc_gga_exc_vxc(const xc_func_type *p, int np, const double *rho, const double *sigma,
		     double *zk, double *vrho, double *vsigma);
void xc_gga_vxc(const xc_func_type *p, int np, const double *rho, const double *sigma,
		 double *vrho, double *vsigma);
void xc_gga_fxc(const xc_func_type *p, int np, const double *rho, const double *sigma,
		 double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_gga_kxc(const xc_func_type *p, int np, const double *rho, const double *sigma,
		 double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);

void xc_gga_lb_modified  (const xc_func_type *p, int np, const double *rho, const double *sigma,
			   double r, double *vrho);

double xc_gga_ak13_get_asymptotic (double homo);

double xc_hyb_exx_coef(const xc_func_type *p);
void  xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta);
void  xc_nlc_coef(const xc_func_type *p, double *nlc_b, double *nlc_C);

/* the meta-GGAs */
void xc_mgga        (const xc_func_type *p, int np,
		      const double *rho, const double *sigma, const double *lapl, const double *tau,
		      double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
		      double *v2rho2, double *v2sigma2, double *v2lapl2, double *v2tau2,
		      double *v2rhosigma, double *v2rholapl, double *v2rhotau,
		      double *v2sigmalapl, double *v2sigmatau, double *v2lapltau);
void xc_mgga_exc    (const xc_func_type *p, int np,
		      const double *rho, const double *sigma, const double *lapl, const double *tau,
		      double *zk);
void xc_mgga_exc_vxc(const xc_func_type *p, int np,
		      const double *rho, const double *sigma, const double *lapl, const double *tau,
		      double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau);
void xc_mgga_vxc    (const xc_func_type *p, int np,
		      const double *rho, const double *sigma, const double *lapl, const double *tau,
		      double *vrho, double *vsigma, double *vlapl, double *vtau);
void xc_mgga_fxc    (const xc_func_type *p, int np,
		      const double *rho, const double *sigma, const double *lapl, const double *tau,
		      double *v2rho2, double *v2sigma2, double *v2lapl2, double *v2tau2,
		      double *v2rhosigma, double *v2rholapl, double *v2rhotau,
		      double *v2sigmalapl, double *v2sigmatau, double *v2lapltau);

/* Functionals that are defined as mixtures of others */
void xc_mix_func
  (const xc_func_type *func, int np,
   const double *rho, const double *sigma, const double *lapl, const double *tau,
   double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
   double *v2rho2, double *v2sigma2, double *v2lapl2, double *v2tau2,
   double *v2rhosigma, double *v2rholapl, double *v2rhotau,
   double *v2sigmalapl, double *v2sigmatau, double *v2lapltau);


#ifdef __cplusplus
}
#endif

#endif
