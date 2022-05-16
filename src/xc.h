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

/* Get the literature reference for libxc */
const char *xc_reference();
/* Get the doi for the literature reference for libxc */
const char *xc_reference_doi();

/* Get the major, minor, and micro version of libxc */
void xc_version(int *major, int *minor, int *micro);
/* Get the version of libxc as a string */
const char *xc_version_string();

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
#define XC_FAMILY_HGGA         32

/* maximum order of derivatives available in libxc */
#define XC_MAXIMUM_ORDER 4

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
#define XC_FLAGS_NEEDS_TAU        (1 << 16) /* 65536 */

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


/* All variables that libxc may input */
#define XC_TOTAL_NUMBER_INPUT_VARIABLES 5

/* spin dimensions of input variables */
typedef union {
  struct {
     int rho, sigma, lapl, tau, exx;
  };
  int fields[XC_TOTAL_NUMBER_INPUT_VARIABLES];
} xc_input_variables_dimensions;

const xc_input_variables_dimensions *input_variables_dimensions_get(int nspin);
  
typedef struct {
  size_t np; /* number of spatial points */
  const xc_input_variables_dimensions *dim; /* spin dimensions of the arrays */

  union {
    struct {
      double *rho;   /* density */
      double *sigma; /* reduced density gradient */
      double *lapl;  /* laplacian of the density */
      double *tau;   /* kinetic energy density */
      double *exx;   /* exchange energy density */
    };
    double *fields[XC_TOTAL_NUMBER_INPUT_VARIABLES];
  };
} xc_input_variables;
  
  
/* All derivatives that libxc may output */
#define XC_TOTAL_NUMBER_OUTPUT_VARIABLES 124

/* spin dimensions of output variables */
typedef union {
  struct {
    int rho, sigma, lapl, tau, exx;       /* spin dimensions of the arrays */
    /* order 0 */
    int zk;
    /* order 1 */
    int vrho, vsigma, vlapl, vtau, vexx;
    /* order 2 */
    int v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2rhoexx;
    int v2sigma2, v2sigmalapl, v2sigmatau, v2sigmaexx;
    int v2lapl2, v2lapltau, v2laplexx;
    int v2tau2, v2tauexx;
    int v2exx2;
    /* order 3 */
    int v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rho2exx;
    int v3rhosigma2, v3rhosigmalapl, v3rhosigmatau, v3rhosigmaexx;
    int v3rholapl2, v3rholapltau, v3rholaplexx;
    int v3rhotau2, v3rhotauexx;
    int v3rhoexx2;
    int v3sigma3, v3sigma2lapl, v3sigma2tau, v3sigma2exx;
    int v3sigmalapl2, v3sigmalapltau, v3sigmalaplexx;
    int v3sigmatau2, v3sigmatauexx;
    int v3sigmaexx2;
    int v3lapl3, v3lapl2tau, v3lapl2exx;
    int v3lapltau2, v3lapltauexx;
    int v3laplexx2;
    int v3tau3, v3tau2exx, v3tauexx2, v3exx3;
    /* order 4 */
    int v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho3exx;
    int v4rho2sigma2, v4rho2sigmalapl, v4rho2sigmatau, v4rho2sigmaexx;
    int v4rho2lapl2, v4rho2lapltau, v4rho2laplexx;
    int v4rho2tau2, v4rho2tauexx;
    int v4rho2exx2;
    int v4rhosigma3, v4rhosigma2lapl, v4rhosigma2tau, v4rhosigma2exx;
    int v4rhosigmalapl2, v4rhosigmalapltau, v4rhosigmalaplexx;
    int v4rhosigmatau2, v4rhosigmatauexx;
    int v4rhosigmaexx2;
    int v4rholapl3, v4rholapl2tau, v4rholapl2exx;
    int v4rholapltau2, v4rholapltauexx;
    int v4rholaplexx2;
    int v4rhotau3, v4rhotau2exx, v4rhoexx3;
    int v4sigma4, v4sigma3lapl, v4sigma3tau, v4sigma3exx;
    int v4sigma2lapl2, v4sigma2lapltau, v4sigma2laplexx;
    int v4sigma2tau2, v4sigma2tauexx;
    int v4sigma2exx2;
    int v4sigmalapl3, v4sigmalapl2tau, v4sigmalapl2exx;
    int v4sigmalapltau2, v4sigmalapltauexx;
    int v4sigmalaplexx2;
    int v4sigmatau3, v4sigmatau2exx, v4sigmatauexx2, v4sigmaexx3;
    int v4lapl4, v4lapl3tau, v4lapl3exx;
    int v4lapl2tau2, v4lapl2tauexx, v4lapl2exx2;
    int v4lapltau3, v4lapltau2exx, v4lapltauexx2, v4laplexx3;
    int v4tau4, v4tau3exx, v4tauexx3, v4exx4;
  };
  int fields[5 + XC_TOTAL_NUMBER_OUTPUT_VARIABLES];
} xc_dimensions;

typedef union { /* this is defined as an union so that we can access the fields sequentially */
  struct {
    /* order 0 (1 var) */
    double *zk;
    
    /* order 1 (5 vars) */
    double *vrho, *vsigma, *vlapl, *vtau, *vexx;
    
    /* order 2 (15 vars) */
    double *v2rho2, *v2rhosigma, *v2rholapl, *v2rhotau, *v2rhoexx;
    double *v2sigma2, *v2sigmalapl, *v2sigmatau, *v2sigmaexx;
    double *v2lapl2, *v2lapltau, *v2laplexx;
    double *v2tau2, *v2tauexx;
    double *v2exx2;
    
    /* order 3 (35 vars) */
    double *v3rho3, *v3rho2sigma, *v3rho2lapl, *v3rho2tau, *v3rho2exx;
    double *v3rhosigma2, *v3rhosigmalapl, *v3rhosigmatau, *v3rhosigmaexx;
    double *v3rholapl2, *v3rholapltau, *v3rholaplexx;
    double *v3rhotau2, *v3rhotauexx;
    double *v3rhoexx2;
    double *v3sigma3, *v3sigma2lapl, *v3sigma2tau, *v3sigma2exx;
    double *v3sigmalapl2, *v3sigmalapltau, *v3sigmalaplexx;
    double *v3sigmatau2, *v3sigmatauexx;
    double *v3sigmaexx2;
    double *v3lapl3, *v3lapl2tau, *v3lapl2exx;
    double *v3lapltau2, *v3lapltauexx;
    double *v3laplexx2;
    double *v3tau3, *v3tau2exx, *v3tauexx2, *v3exx3;
    
    /* order 4 (68 vars) */
    double *v4rho4, *v4rho3sigma, *v4rho3lapl, *v4rho3tau, *v4rho3exx;
    double *v4rho2sigma2, *v4rho2sigmalapl, *v4rho2sigmatau, *v4rho2sigmaexx;
    double *v4rho2lapl2, *v4rho2lapltau, *v4rho2laplexx;
    double *v4rho2tau2, *v4rho2tauexx;
    double *v4rho2exx2;
    double *v4rhosigma3, *v4rhosigma2lapl, *v4rhosigma2tau, *v4rhosigma2exx;
    double *v4rhosigmalapl2, *v4rhosigmalapltau, *v4rhosigmalaplexx;
    double *v4rhosigmatau2, *v4rhosigmatauexx;
    double *v4rhosigmaexx2;
    double *v4rholapl3, *v4rholapl2tau, *v4rholapl2exx;
    double *v4rholapltau2, *v4rholapltauexx;
    double *v4rholaplexx2;
    double *v4rhotau3, *v4rhotau2exx, *v4rhoexx3;
    double *v4sigma4, *v4sigma3lapl, *v4sigma3tau, *v4sigma3exx;
    double *v4sigma2lapl2, *v4sigma2lapltau, *v4sigma2laplexx;
    double *v4sigma2tau2, *v4sigma2tauexx;
    double *v4sigma2exx2;
    double *v4sigmalapl3, *v4sigmalapl2tau, *v4sigmalapl2exx;
    double *v4sigmalapltau2, *v4sigmalapltauexx;
    double *v4sigmalaplexx2;
    double *v4sigmatau3, *v4sigmatau2exx, *v4sigmatauexx2, *v4sigmaexx3;
    double *v4lapl4, *v4lapl3tau, *v4lapl3exx;
    double *v4lapl2tau2, *v4lapl2tauexx, *v4lapl2exx2;
    double *v4lapltau3, *v4lapltau2exx, *v4lapltauexx2, *v4laplexx3;
    double *v4tau4, *v4tau3exx, *v4tauexx3, *v4exx4;
  };
  double *fields[XC_TOTAL_NUMBER_OUTPUT_VARIABLES];
} xc_output_variables;

/* from io_variables.c */
extern const char *xc_input_variables_name[];     /* mapping input variable -> name */
extern const int xc_input_variables_family_key[]; /* mapping input variable -> family */
extern const int xc_input_variables_flags_key[];  /* mapping input variable -> flags */

xc_input_variables *xc_input_variables_allocate(double np, int family, int flags, int nspin);
int xc_input_variables_sanity_check(const xc_input_variables *out, int family, int flags);
void xc_input_variables_initialize(xc_input_variables *out);
void xc_input_variables_deallocate(xc_input_variables *out);
  
extern const char *xc_output_variables_name[];     /* mapping output variable -> name */
extern const int xc_output_variables_order_key[];  /* mapping output variable -> order of derivative */
extern const int xc_output_variables_family_key[]; /* mapping output variable -> family */
extern const int xc_output_variables_flags_key[];  /* mapping output variable -> flags */

xc_output_variables *xc_output_variables_allocate(double np, const int *orders, int family, int flags, int nspin);
int xc_output_variables_sanity_check(const xc_output_variables *out, const int *orders, int family, int flags);
void xc_output_variables_initialize(xc_output_variables *out, int np, int nspin);
void xc_output_variables_deallocate(xc_output_variables *out);

/* type of the lda function */
typedef void (*xc_functionals_work)(const struct xc_func_type *p,
     const xc_input_variables *in, xc_output_variables *out);
  
typedef struct {
  const xc_functionals_work unpol[5], pol[5];
} xc_functionals_work_variants;
  
/* Structure that contains information on the functional */
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
  const xc_functionals_work_variants  *work;
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

  const xc_dimensions *dim;    /* the dimensions of all input and output arrays */

  /* This is where the values of the external parameters are stored */
  double *ext_params;
  /* This is a placeholder for structs of parameters that are used in the Maple generated sources */
  void *params;

  double dens_threshold;       /* functional is put to zero for spin-densities smaller than this */
  double zeta_threshold;       /* idem for the absolute value of zeta */
  double sigma_threshold;
  double tau_threshold;
};

typedef struct xc_func_type xc_func_type;


/** Get a functional's id number from its name  */
int   xc_functional_get_number(const char *name);
/** Get a functional's name from its id number  */
char *xc_functional_get_name(int number);
/** Get a functional's family and the number within the family from the id number */
int   xc_family_from_id(int id, int *family, int *number);

/** The number of functionals implemented in this version of libxc */
int   xc_number_of_functionals();
/** The maximum name length of any functional */
int   xc_maximum_name_length();
/** Returns the available functional number sorted by id */
void  xc_available_functional_numbers(int *list);
/** Returns the available functional number sorted by the functionals'
    names; this function is a helper for the Python frontend. */
void  xc_available_functional_numbers_by_name(int *list);
/** Fills the list with the names of the available functionals,
    ordered by name. The list array should be [Nfuncs][maxlen+1]. */
void  xc_available_functional_names(char **list);

/** Dynamically allocates a libxc functional; which will also need to be initialized. */
xc_func_type *xc_func_alloc();
/** Initializes a functional by id with nspin spin channels */
int   xc_func_init(xc_func_type *p, int functional, int nspin);
/** Destructor for an initialized functional */
void  xc_func_end(xc_func_type *p);
/** Frees a dynamically allocated functional */
void  xc_func_free(xc_func_type *p);
/** Get information on a functional */
const xc_func_info_type *xc_func_get_info(const xc_func_type *p);

/** Sets the density threshold for a functional */
void  xc_func_set_dens_threshold(xc_func_type *p, double t_dens);
/** Sets the spin polarization threshold for a functional */
void  xc_func_set_zeta_threshold(xc_func_type *p, double t_zeta);
/** Sets the reduced gradient threshold for a functional */
void  xc_func_set_sigma_threshold(xc_func_type *p, double t_sigma);
/** Sets the kinetic energy density threshold for a functional */
void  xc_func_set_tau_threshold(xc_func_type *p, double t_tau);

/** Sets all external parameters for a functional */
void  xc_func_set_ext_params(xc_func_type *p, const double *ext_params);
/** Gets all external parameters for a functional. Array needs to be preallocated  */
void  xc_func_get_ext_params(const xc_func_type *p, double *ext_params);
/** Sets an external parameter by name for a functional */
void  xc_func_set_ext_params_name(xc_func_type *p, const char *name, double par);
/** Gets an external parameter by name for a functional */
double xc_func_get_ext_params_name(const xc_func_type *p, const char *name);
/** Gets an external parameter by index */
double xc_func_get_ext_params_value(const xc_func_type *p, int number);

#include "xc_funcs.h"
#include "xc_funcs_removed.h"

/* New API */
void xc_evaluate_func(const xc_func_type *p, int order,
                      const xc_input_variables *in, xc_output_variables *out);

/* This is the old, Fortran friendly interface */

/*
  the LDAs
*/
void xc_lda(const xc_func_type *p, size_t np, double *rho,
       double *zk, double *vrho, double *v2rho2, double *v3rho3, double *v4rho4);
void xc_lda_exc(const xc_func_type *p, size_t np, double *rho,
       double *zk);
void xc_lda_exc_vxc(const xc_func_type *p, size_t np, double *rho,
       double *zk, double *vrho);
void xc_lda_exc_vxc(const xc_func_type *p, size_t np, double *rho,
       double *zk, double *vrho);
void xc_lda_exc_vxc_fxc(const xc_func_type *p, size_t np, double *rho,
       double *zk, double *vrho, double *v2rho2);
void xc_lda_vxc_fxc(const xc_func_type *p, size_t np, double *rho,
       double *vrho, double *v2rho2);
void xc_lda_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho,
       double *zk, double *vrho, double *v2rho2, double *v3rho3);
void xc_lda_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho,
       double *vrho, double *v2rho2, double *v3rho3);
void xc_lda_vxc(const xc_func_type *p, size_t np, double *rho,
       double *vrho);
void xc_lda_fxc(const xc_func_type *p, size_t np, double *rho,
       double *v2rho2);
void xc_lda_kxc(const xc_func_type *p, size_t np, double *rho,
       double *v3rho3);
void xc_lda_lxc(const xc_func_type *p, size_t np, double *rho,
       double *v4rho4);
  
/*
  the GGAs
*/
void xc_gga(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2,
       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3,
       double *v4rho4, double *v4rho3sigma, double *v4rho2sigma2, double *v4rhosigma3, double *v4sigma4);
void xc_gga_exc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk);
void xc_gga_exc_vxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk,
       double *vrho, double *vsigma);
void xc_gga_exc_vxc_fxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_gga_vxc_fxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_gga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *zk,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2,
       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_gga_vxc_fxc_kxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *vrho, double *vsigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2,
       double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_gga_vxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
		   double *vrho, double *vsigma);
void xc_gga_fxc(const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *v2rho2, double *v2rhosigma, double *v2sigma2);
void xc_gga_kxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
  		 double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3);
void xc_gga_lxc (const xc_func_type *p, size_t np, double *rho, double *sigma,
       double *v4rho4,  double *v4rho3sigma,  double *v4rho2sigma2,  double *v4rhosigma3,
       double *v4sigma4);

/*
  the mGGAs
*/
void xc_mgga(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau, double *v2sigma2,
       double *v2sigmalapl, double *v2sigmatau, double *v2lapl2, double *v2lapltau, double *v2tau2,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau, double *v3rhosigma2,
       double *v3rhosigmalapl, double *v3rhosigmatau, double *v3rholapl2, double *v3rholapltau,
       double *v3rhotau2, double *v3sigma3, double *v3sigma2lapl, double *v3sigma2tau,
       double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmatau2, double *v3lapl3,
       double *v3lapl2tau, double *v3lapltau2, double *v3tau3,
       double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
       double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
       double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
       double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
       double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
       double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
       double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
       double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
       double *v4lapl2tau2, double *v4lapltau3, double *v4tau4);
void xc_mgga_exc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *zk);
void xc_mgga_exc_vxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau);
void xc_mgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
       double *v2lapltau, double *v2tau2);
void xc_mgga_vxc_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *vrho, double *vsigma, double *vlapl, double *vtau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
       double *v2lapltau, double *v2tau2);
void xc_mgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
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
void xc_mgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
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
void xc_mgga_vxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *vrho, double *vsigma, double *vlapl, double *vtau);
void xc_mgga_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2sigma2, double *v2sigmalapl, double *v2sigmatau, double *v2lapl2,
       double *v2lapltau, double *v2tau2);
void xc_mgga_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
       double *v3rhosigma2, double *v3rhosigmalapl, double *v3rhosigmatau,
       double *v3rholapl2, double *v3rholapltau, double *v3rhotau2, double *v3sigma3,
       double *v3sigma2lapl, double *v3sigma2tau, double *v3sigmalapl2, double *v3sigmalapltau,
       double *v3sigmatau2, double *v3lapl3, double *v3lapl2tau, double *v3lapltau2,
       double *v3tau3);
void xc_mgga_lxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau,
       double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau, double *v4rho2sigma2,
       double *v4rho2sigmalapl, double *v4rho2sigmatau, double *v4rho2lapl2, double *v4rho2lapltau,
       double *v4rho2tau2, double *v4rhosigma3, double *v4rhosigma2lapl, double *v4rhosigma2tau,
       double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmatau2,
       double *v4rholapl3, double *v4rholapl2tau, double *v4rholapltau2, double *v4rhotau3,
       double *v4sigma4, double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma2lapl2,
       double *v4sigma2lapltau, double *v4sigma2tau2, double *v4sigmalapl3, double *v4sigmalapl2tau,
       double *v4sigmalapltau2, double *v4sigmatau3, double *v4lapl4, double *v4lapl3tau,
       double *v4lapl2tau2, double *v4lapltau3, double *v4tau4);

  /* 
     the HGGAs 
  */
void xc_hgga(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
       double *v3rho2exx, double *v3rhosigma2, double *v3rhosigmalapl,
       double *v3rhosigmatau, double *v3rhosigmaexx, double *v3rholapl2,
       double *v3rholapltau, double *v3rholaplexx, double *v3rhotau2,
       double *v3rhotauexx, double *v3rhoexx2, double *v3sigma3, double *v3sigma2lapl,
       double *v3sigma2tau, double *v3sigma2exx, double *v3sigmalapl2,
       double *v3sigmalapltau, double *v3sigmalaplexx, double *v3sigmatau2,
       double *v3sigmatauexx, double *v3sigmaexx2, double *v3lapl3,
       double *v3lapl2tau, double *v3lapl2exx, double *v3lapltau2,
       double *v3lapltauexx, double *v3laplexx2, double *v3tau3, double *v3tau2exx,
       double *v3tauexx2, double *v3exx3,
       double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau,
       double *v4rho3exx, double *v4rho2sigma2, double *v4rho2sigmalapl,
       double *v4rho2sigmatau, double *v4rho2sigmaexx, double *v4rho2lapl2,
       double *v4rho2lapltau, double *v4rho2laplexx, double *v4rho2tau2,
       double *v4rho2tauexx, double *v4rho2exx2, double *v4rhosigma3,
       double *v4rhosigma2lapl, double *v4rhosigma2tau, double *v4rhosigma2exx,
       double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmalaplexx,
       double *v4rhosigmatau2, double *v4rhosigmatauexx, double *v4rhosigmaexx2,
       double *v4rholapl3, double *v4rholapl2tau, double *v4rholapl2exx,
       double *v4rholapltau2, double *v4rholapltauexx, double *v4rholaplexx2,
       double *v4rhotau3, double *v4rhotau2exx, double *v4rhoexx3, double *v4sigma4,
       double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma3exx,
       double *v4sigma2lapl2, double *v4sigma2lapltau, double *v4sigma2laplexx,
       double *v4sigma2tau2, double *v4sigma2tauexx, double *v4sigma2exx2,
       double *v4sigmalapl3, double *v4sigmalapl2tau, double *v4sigmalapl2exx,
       double *v4sigmalapltau2, double *v4sigmalapltauexx, double *v4sigmalaplexx2,
       double *v4sigmatau3, double *v4sigmatau2exx, double *v4sigmatauexx2,
       double *v4sigmaexx3, double *v4lapl4, double *v4lapl3tau, double *v4lapl3exx,
       double *v4lapl2tau2, double *v4lapl2tauexx, double *v4lapl2exx2,
       double *v4lapltau3, double *v4lapltau2exx, double *v4lapltauexx2,
       double *v4laplexx3, double *v4tau4, double *v4tau3exx, double *v4tauexx3,
       double *v4exx4);
void xc_hgga_exc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *zk);
void xc_hgga_exc_vxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx);
void xc_hgga_exc_vxc_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2);
void xc_hgga_vxc_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2);
void xc_hgga_exc_vxc_fxc_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *zk,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
       double *v3rho2exx, double *v3rhosigma2, double *v3rhosigmalapl,
       double *v3rhosigmatau, double *v3rhosigmaexx, double *v3rholapl2,
       double *v3rholapltau, double *v3rholaplexx, double *v3rhotau2,
       double *v3rhotauexx, double *v3rhoexx2, double *v3sigma3, double *v3sigma2lapl,
       double *v3sigma2tau, double *v3sigma2exx, double *v3sigmalapl2,
       double *v3sigmalapltau, double *v3sigmalaplexx, double *v3sigmatau2,
       double *v3sigmatauexx, double *v3sigmaexx2, double *v3lapl3,
       double *v3lapl2tau, double *v3lapl2exx, double *v3lapltau2,
       double *v3lapltauexx, double *v3laplexx2, double *v3tau3, double *v3tau2exx,
       double *v3tauexx2, double *v3exx3);
void xc_hgga_vxc_fxc_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
       double *v3rho2exx, double *v3rhosigma2, double *v3rhosigmalapl,
       double *v3rhosigmatau, double *v3rhosigmaexx, double *v3rholapl2,
       double *v3rholapltau, double *v3rholaplexx, double *v3rhotau2,
       double *v3rhotauexx, double *v3rhoexx2, double *v3sigma3, double *v3sigma2lapl,
       double *v3sigma2tau, double *v3sigma2exx, double *v3sigmalapl2,
       double *v3sigmalapltau, double *v3sigmalaplexx, double *v3sigmatau2,
       double *v3sigmatauexx, double *v3sigmaexx2, double *v3lapl3,
       double *v3lapl2tau, double *v3lapl2exx, double *v3lapltau2,
       double *v3lapltauexx, double *v3laplexx2, double *v3tau3, double *v3tau2exx,
       double *v3tauexx2, double *v3exx3);
void xc_hgga_vxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *vrho, double *vsigma, double *vlapl, double *vtau, double *vexx);
void xc_hgga_fxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *v2rho2, double *v2rhosigma, double *v2rholapl, double *v2rhotau,
       double *v2rhoexx, double *v2sigma2, double *v2sigmalapl, double *v2sigmatau,
       double *v2sigmaexx, double *v2lapl2, double *v2lapltau, double *v2laplexx,
       double *v2tau2, double *v2tauexx, double *v2exx2);
void xc_hgga_mgga_kxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *v3rho3, double *v3rho2sigma, double *v3rho2lapl, double *v3rho2tau,
       double *v3rho2exx, double *v3rhosigma2, double *v3rhosigmalapl,
       double *v3rhosigmatau, double *v3rhosigmaexx, double *v3rholapl2,
       double *v3rholapltau, double *v3rholaplexx, double *v3rhotau2,
       double *v3rhotauexx, double *v3rhoexx2, double *v3sigma3,
       double *v3sigma2lapl, double *v3sigma2tau, double *v3sigma2exx,
       double *v3sigmalapl2, double *v3sigmalapltau, double *v3sigmalaplexx,
       double *v3sigmatau2, double *v3sigmatauexx, double *v3sigmaexx2,
       double *v3lapl3, double *v3lapl2tau, double *v3lapl2exx, double *v3lapltau2,
       double *v3lapltauexx, double *v3laplexx2, double *v3tau3, double *v3tau2exx,
       double *v3tauexx2, double *v3exx3);
void xc_hgga_lxc(const xc_func_type *p, size_t np,
       double *rho, double *sigma, double *lapl, double *tau, double *exx,
       double *v4rho4, double *v4rho3sigma, double *v4rho3lapl, double *v4rho3tau,
       double *v4rho3exx, double *v4rho2sigma2, double *v4rho2sigmalapl,
       double *v4rho2sigmatau, double *v4rho2sigmaexx, double *v4rho2lapl2,
       double *v4rho2lapltau, double *v4rho2laplexx, double *v4rho2tau2,
       double *v4rho2tauexx, double *v4rho2exx2, double *v4rhosigma3,
       double *v4rhosigma2lapl, double *v4rhosigma2tau, double *v4rhosigma2exx,
       double *v4rhosigmalapl2, double *v4rhosigmalapltau, double *v4rhosigmalaplexx,
       double *v4rhosigmatau2, double *v4rhosigmatauexx, double *v4rhosigmaexx2,
       double *v4rholapl3, double *v4rholapl2tau, double *v4rholapl2exx,
       double *v4rholapltau2, double *v4rholapltauexx, double *v4rholaplexx2,
       double *v4rhotau3, double *v4rhotau2exx, double *v4rhoexx3, double *v4sigma4,
       double *v4sigma3lapl, double *v4sigma3tau, double *v4sigma3exx,
       double *v4sigma2lapl2, double *v4sigma2lapltau, double *v4sigma2laplexx,
       double *v4sigma2tau2, double *v4sigma2tauexx, double *v4sigma2exx2,
       double *v4sigmalapl3, double *v4sigmalapl2tau, double *v4sigmalapl2exx,
       double *v4sigmalapltau2, double *v4sigmalapltauexx, double *v4sigmalaplexx2,
       double *v4sigmatau3, double *v4sigmatau2exx, double *v4sigmatauexx2,
       double *v4sigmaexx3, double *v4lapl4, double *v4lapl3tau, double *v4lapl3exx,
       double *v4lapl2tau2, double *v4lapl2tauexx, double *v4lapl2exx2,
       double *v4lapltau3, double *v4lapltau2exx, double *v4lapltauexx2,
       double *v4laplexx3, double *v4tau4, double *v4tau3exx, double *v4tauexx3,
       double *v4exx4);

  
/* Calculate asymptotic value of the AK13 potential */
double xc_gga_ak13_get_asymptotic (double homo);
/* Calculate asymptotic value of the AK13 potential with customized parameter values */
double xc_gga_ak13_pars_get_asymptotic (double homo, const double *ext_params);

/* Returns the hybrid type of a functional */
int xc_hyb_type(const xc_func_type *p);
/* Returns fraction of Hartree-Fock exchange in a global hybrid functional */
double xc_hyb_exx_coef(const xc_func_type *p);
/* Returns fraction of Hartee-Fock exchange and short-range exchange in a range-separated hybrid functional  */
void xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta);
/* Returns the b and C coefficients for a non-local VV10 correlation kernel */
void xc_nlc_coef(const xc_func_type *p, double *nlc_b, double *nlc_C);

/* If this is a mixed functional, returns the number of auxiliary functions. Otherwise returns zero. */
int xc_num_aux_funcs(const xc_func_type *p);
/* Gets the IDs of the auxiliary functions */
void xc_aux_func_ids(const xc_func_type *p, int *ids);
/* Gets the weights of the auxiliary functions */
void xc_aux_func_weights(const xc_func_type *p, double *weights);

#ifdef __cplusplus
}
#endif

#endif
