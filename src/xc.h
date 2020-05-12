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

#include <xc_version.h>
#include <stddef.h>

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
#define XC_FLAGS_VV10             (1 << 10) /*  1024 */
#define XC_FLAGS_STABLE           (1 << 13) /*  8192 */
/* functionals marked with the development flag may have significant problems in the implementation */
#define XC_FLAGS_DEVELOPMENT      (1 << 14) /* 16384 */
#define XC_FLAGS_NEEDS_LAPLACIAN  (1 << 15) /* 32768 */

/* This is the case for most functionals in libxc */
#define XC_FLAGS_HAVE_ALL         (XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | \
                                   XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC | \
                                   XC_FLAGS_HAVE_LXC)

/* This magic value means use default parameter */
#define XC_EXT_PARAMS_DEFAULT   -999998888

/* Different flavors of many-body terms used in hybrids
   The Fock term to be added to the Hamiltonian reads

     F = -1/2 <i j| f(r_12)/r_12 |j i>

   where the function f(r) is

   *) XC_HYB_FOCK           f(r) = coeff
   *) XC_HYB_ERF_SR         f(r) = coeff * (1 - erf(omega r))
   *) XC_HYB_YUKAWA_SR      f(r) = coeff * exp(-omega r)
   *) XC_HYB_GAUSSIAN_SR    f(r) = coeff * 2*omega/sqrt(pi) * exp(-omega^2 r^2)
*/
#define XC_HYB_NONE             0
#define XC_HYB_FOCK             1  /* Normal hybrid */
#define XC_HYB_PT2              2  /* Used for double hybrids */
#define XC_HYB_ERF_SR           4  /* Short range of range separated - erf version */
#define XC_HYB_YUKAWA_SR        8  /* Short range of range separated - Yakawa version */
#define XC_HYB_GAUSSIAN_SR     16  /* Short range of range separated - Gaussian version */

/* Different types of hybrid functionals. */
#define XC_HYB_SEMILOCAL        0  /* Standard semi-local functional (not a hybrid) */
#define XC_HYB_HYBRID           1  /* Standard hybrid functional */
#define XC_HYB_CAM              2  /* Coulomb attenuated hybrid */
#define XC_HYB_CAMY             3  /* Coulomb attenuated hybrid with a Yukawa screening */
#define XC_HYB_CAMG             4  /* Coulomb attenuated hybrid with a Gaussian screening */
#define XC_HYB_DOUBLE_HYBRID    5  /* Double hybrid */
#define XC_HYB_MIXTURE      32768  /* More complicated mixture (have to check individual terms) */
  
#define XC_MAX_REFERENCES       5

/* This are the derivatives that a functional returns */
#define XC_NOARG
#define XC_COMMA ,
#define LDA_OUT_PARAMS_NO_EXC(P1_, P2_) \
  P1_ P2_ ## vrho   P1_ P2_ ## v2rho2 \
  P1_ P2_ ## v3rho3 P1_ P2_ ## v4rho4

#define GGA_OUT_PARAMS_NO_EXC(P1_, P2_) \
  P1_ P2_ ## vrho         P1_ P2_ ## vsigma       \
  P1_ P2_ ## v2rho2       P1_ P2_ ## v2rhosigma   \
  P1_ P2_ ## v2sigma2                             \
  P1_ P2_ ## v3rho3       P1_ P2_ ## v3rho2sigma  \
  P1_ P2_ ## v3rhosigma2  P1_ P2_ ## v3sigma3     \
  P1_ P2_ ## v4rho4       P1_ P2_ ## v4rho3sigma  \
  P1_ P2_ ## v4rho2sigma2 P1_ P2_ ## v4rhosigma3  \
  P1_ P2_ ## v4sigma4

/* This are the derivatives of a mgga
       1st order:  4
       2nd order: 10
       3rd order: 20
       4th order: 35
 */
#define MGGA_OUT_PARAMS_NO_EXC(P1_, P2_) \
  P1_ P2_ ## vrho              P1_ P2_ ## vsigma          \
  P1_ P2_ ## vlapl             P1_ P2_ ## vtau            \
  P1_ P2_ ## v2rho2            P1_ P2_ ## v2rhosigma      \
  P1_ P2_ ## v2rholapl         P1_ P2_ ## v2rhotau        \
  P1_ P2_ ## v2sigma2          P1_ P2_ ## v2sigmalapl     \
  P1_ P2_ ## v2sigmatau        P1_ P2_ ## v2lapl2         \
  P1_ P2_ ## v2lapltau         P1_ P2_ ## v2tau2          \
  P1_ P2_ ## v3rho3            P1_ P2_ ## v3rho2sigma     \
  P1_ P2_ ## v3rho2lapl        P1_ P2_ ## v3rho2tau       \
  P1_ P2_ ## v3rhosigma2       P1_ P2_ ## v3rhosigmalapl  \
  P1_ P2_ ## v3rhosigmatau     P1_ P2_ ## v3rholapl2      \
  P1_ P2_ ## v3rholapltau      P1_ P2_ ## v3rhotau2       \
  P1_ P2_ ## v3sigma3          P1_ P2_ ## v3sigma2lapl    \
  P1_ P2_ ## v3sigma2tau       P1_ P2_ ## v3sigmalapl2    \
  P1_ P2_ ## v3sigmalapltau    P1_ P2_ ## v3sigmatau2     \
  P1_ P2_ ## v3lapl3           P1_ P2_ ## v3lapl2tau      \
  P1_ P2_ ## v3lapltau2        P1_ P2_ ## v3tau3          \
  P1_ P2_ ## v4rho4            P1_ P2_ ## v4rho3sigma     \
  P1_ P2_ ## v4rho3lapl        P1_ P2_ ## v4rho3tau       \
  P1_ P2_ ## v4rho2sigma2      P1_ P2_ ## v4rho2sigmalapl \
  P1_ P2_ ## v4rho2sigmatau    P1_ P2_ ## v4rho2lapl2     \
  P1_ P2_ ## v4rho2lapltau     P1_ P2_ ## v4rho2tau2      \
  P1_ P2_ ## v4rhosigma3       P1_ P2_ ## v4rhosigma2lapl \
  P1_ P2_ ## v4rhosigma2tau    P1_ P2_ ## v4rhosigmalapl2 \
  P1_ P2_ ## v4rhosigmalapltau P1_ P2_ ## v4rhosigmatau2  \
  P1_ P2_ ## v4rholapl3        P1_ P2_ ## v4rholapl2tau   \
  P1_ P2_ ## v4rholapltau2     P1_ P2_ ## v4rhotau3       \
  P1_ P2_ ## v4sigma4          P1_ P2_ ## v4sigma3lapl    \
  P1_ P2_ ## v4sigma3tau       P1_ P2_ ## v4sigma2lapl2   \
  P1_ P2_ ## v4sigma2lapltau   P1_ P2_ ## v4sigma2tau2    \
  P1_ P2_ ## v4sigmalapl3      P1_ P2_ ## v4sigmalapl2tau \
  P1_ P2_ ## v4sigmalapltau2   P1_ P2_ ## v4sigmatau3     \
  P1_ P2_ ## v4lapl4           P1_ P2_ ## v4lapl3tau      \
  P1_ P2_ ## v4lapl2tau2       P1_ P2_ ## v4lapltau3      \
  P1_ P2_ ## v4tau4


void xc_version(int *major, int *minor, int *micro);
const char *xc_version_string();

struct xc_func_type;

typedef struct{
  const char *ref, *doi, *bibtex;
} func_reference_type;

char const *xc_func_reference_get_ref(const func_reference_type *reference);
char const *xc_func_reference_get_doi(const func_reference_type *reference);
char const *xc_func_reference_get_bibtex(const func_reference_type *reference);

  
typedef struct{
  int n; /* Number of parameters */

  const char **names; /* ATTENTION: if name starts with a _ it is an *internal* parameter,
                        changing the value effectively changes the functional! */
  const char **descriptions; /* long description of the parameters */
  const double *values; /* default values of the parameters */

  void (*set)(struct xc_func_type *p, const double *ext_params);
} func_params_type;

  
typedef struct{
  int   number;   /* identifier number */
  int   kind;     /* XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC */

  const char *name;     /* name of the functional, e.g. "PBE" */
  int   family;   /* type of the functional, e.g. XC_FAMILY_GGA */
  func_reference_type *refs[XC_MAX_REFERENCES];  /* index of the references */

  int   flags;    /* see above for a list of possible flags */

  double dens_threshold;

  /* this allows to have external parameters in the functional */
  func_params_type ext_params;

  void (*init)(struct xc_func_type *p);
  void (*end) (struct xc_func_type *p);
  void (*lda) (const struct xc_func_type *p, size_t np,
               const double *rho,
               double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
  void (*gga) (const struct xc_func_type *p, size_t np,
               const double *rho, const double *sigma,
               double *zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
  void (*mgga)(const struct xc_func_type *p, size_t np,
               const double *rho, const double *sigma, const double *lapl_rho, const double *tau,
               double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
} xc_func_info_type;

  
/* for API compability with older versions of libxc */
#define XC(func) xc_ ## func

  
int xc_func_info_get_number(const xc_func_info_type *info);
int xc_func_info_get_kind(const xc_func_info_type *info);
char const *xc_func_info_get_name(const xc_func_info_type *info);
int xc_func_info_get_family(const xc_func_info_type *info);
int xc_func_info_get_flags(const xc_func_info_type *info);
const func_reference_type *xc_func_info_get_references(const xc_func_info_type *info, int number);

  
int xc_func_info_get_n_ext_params(const xc_func_info_type *info);
char const *xc_func_info_get_ext_params_name(const xc_func_info_type *p, int number);
char const *xc_func_info_get_ext_params_description(const xc_func_info_type *info, int number);
double xc_func_info_get_ext_params_default_value(const xc_func_info_type *info, int number);

  
struct xc_dimensions{
  int rho, sigma, lapl, tau;       /* spin dimensions of the arrays */
  int zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA, );
};

typedef struct xc_dimensions xc_dimensions;

  
struct xc_func_type{
  const xc_func_info_type *info;       /* all the information concerning this functional */
  int nspin;                           /* XC_UNPOLARIZED or XC_POLARIZED  */

  int n_func_aux;                      /* how many auxiliary functions we need */
  struct xc_func_type **func_aux;      /* most GGAs are based on a LDA or other GGAs  */
  double *mix_coef;                    /* coefficients for the mixing */

  /**
     Parameters for range-separated hybrids
     hyb_type[i]:  XC_HYB_NONE, XC_HYB_FOCK, XC_HYB_ERF_SR, etc.
     hyb_omega[i]: the range separation constant
     hyb_coeff[i]: fraction of exchange, used both for
                usual hybrids as well as range-separated ones

     N.B. Different conventions for alpha and beta can be found in
     literature. In the convention used in libxc, at short range the
     fraction of exact exchange is cam_alpha+cam_beta, while at long
     range it is cam_alpha.
  */
  int hyb_number_terms, *hyb_type;
  double *hyb_coeff, *hyb_omega;

  double nlc_b;                /* Non-local correlation, b parameter */
  double nlc_C;                /* Non-local correlation, C parameter */

  xc_dimensions dim;           /* the dimensions of all input and output arrays */

  void *params;                /* this allows us to fix parameters in the functional */
  double dens_threshold;
};

typedef struct xc_func_type xc_func_type;

  
/* functionals */
int   xc_functional_get_number(const char *name);
char *xc_functional_get_name(int number);
int   xc_family_from_id(int id, int *family, int *number);
int   xc_number_of_functionals();
int   xc_maximum_name_length();
void  xc_available_functional_numbers(int *list);
void  xc_available_functional_names(char **list);

xc_func_type *xc_func_alloc();
int   xc_func_init(xc_func_type *p, int functional, int nspin);
void  xc_func_end(xc_func_type *p);
void  xc_func_free(xc_func_type *p);
const xc_func_info_type *xc_func_get_info(const xc_func_type *p);
void  xc_func_set_dens_threshold(xc_func_type *p, double dens_threshold);
void  xc_func_set_ext_params(xc_func_type *p, const double *ext_params);
void  xc_func_set_ext_params_name(xc_func_type *p, const char *name, double par);

  
#include "xc_funcs.h"
#include "xc_funcs_removed.h"

  
void xc_lda (const xc_func_type *p, size_t np, const double *rho,
             double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
void xc_gga (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
             double *zk GGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));
void xc_mgga(const xc_func_type *p, size_t np,
             const double *rho, const double *sigma, const double *lapl_rho, const double *tau,
             double *zk MGGA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ));

void xc_lda_exc (const xc_func_type *p, size_t np, const double *rho, double *zk);
void xc_gga_exc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		 double *zk);
void xc_mgga_exc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *zk);

#ifndef XC_DONT_COMPILE_VXC
void xc_lda_exc_vxc (const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho);
void xc_gga_exc_vxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		 double *zk, double *vrho, double *vsigma);
void xc_mgga_exc_vxc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau);

void xc_lda_vxc (const xc_func_type *p, size_t np, const double *rho, double *vrho);
void xc_gga_vxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		 double *vrho, double *vsigma);
void xc_mgga_vxc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *vrho, double *vsigma, double *vlapl, double *vtau);

#ifndef XC_DONT_COMPILE_FXC
void xc_lda_exc_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2);
void xc_gga_exc_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                         double *zk, double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_mgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2);

void xc_lda_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2);
void xc_gga_vxc_fxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                         double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_mgga_vxc_fxc(const xc_func_type *p, size_t np,
                         const double *rho, const double *sigma, const double *lapl, const double *tau,
                         double *vrho, double *vsigma, double *vlapl, double *vtau,
                         double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                         double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                         double *v2lapltau, double *v2tau2);

void xc_lda_fxc (const xc_func_type *p, size_t np, const double *rho, double *v2rho2);
void xc_gga_fxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		 double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_mgga_fxc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
     double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
     double *v2lapltau, double *v2tau2);

#ifndef XC_DONT_COMPILE_KXC
void xc_lda_exc_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3);
void xc_gga_exc_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                             double *zk, double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                             double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_mgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *zk, double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3);

void xc_lda_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, double *vrho, double *v2rho2, double *v3rho3);
void xc_gga_vxc_fxc_kxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
                             double *vrho, double *vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2,
                             double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_mgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
                             const double *rho, const double *sigma, const double *lapl, const double *tau,
                             double *vrho, double *vsigma, double *vlapl, double *vtau,
                             double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
                             double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
                             double *v2lapltau, double *v2tau2,
                             double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
                             double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
                             double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
                             double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
                             double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
                             double *v3tau3);

void xc_lda_kxc (const xc_func_type *p, size_t np, const double *rho, double *v3rho3);
void xc_gga_kxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
		 double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_mgga_kxc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
     double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
     double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
     double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
     double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
     double *v3tau3);

#ifndef XC_DONT_COMPILE_LXC
void xc_lda_lxc (const xc_func_type *p, size_t np, const double *rho, double *v4rho4);
void xc_gga_lxc (const xc_func_type *p, size_t np, const double *rho, const double *sigma,
     double *v4rho4,  double *v4rho3sigma,  double *v4rho2sigma2,  double *v4rhosigma3,
     double *v4sigma4);
void xc_mgga_lxc(const xc_func_type *p, size_t np,
     const double *rho, const double *sigma, const double *lapl, const double *tau,
     double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
     double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
     double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
     double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
     double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
     double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
     double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
     double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
     double *v4lapl2tau2, double *v4lapltau3, double *v4tau4);
#endif
#endif
#endif
#endif

/* Calculate asymptotic value of the AK13 potential */
double xc_gga_ak13_get_asymptotic (double homo);
/* Calculate asymptotic value of the AK13 potential with customized parameter values */
double xc_gga_ak13_pars_get_asymptotic (double homo, const double *ext_params);


/* the hybrids and van der Waals functions */
int xc_hyb_type(const xc_func_type *p);
double xc_hyb_exx_coef(const xc_func_type *p);
void xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta);
void xc_nlc_coef(const xc_func_type *p, double *nlc_b, double *nlc_C);

#ifdef __cplusplus
}
#endif

#endif
