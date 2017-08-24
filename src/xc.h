/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
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
#include "xc_config.h"
  
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
  
void XC(version)(int *major, int *minor, int *micro);
const char *XC(version_string)();
    
struct XC(func_type);

typedef struct{
  char *ref, *doi, *bibtex;
} func_reference_type;

char const *XC(func_reference_get_ref)(const func_reference_type *reference);
char const *XC(func_reference_get_doi)(const func_reference_type *reference);
char const *XC(func_reference_get_bibtex)(const func_reference_type *reference);

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

  FLOAT min_dens, min_grad, min_tau;

  /* this allows to have external parameters in the functional */
  int n_ext_params;
  const func_params_type *ext_params;
  void (*set_ext_params)(struct XC(func_type) *p, const double *ext_params);

  void (*init)(struct XC(func_type) *p);
  void (*end) (struct XC(func_type) *p);
  void (*lda) (const struct XC(func_type) *p, int np, 
	       const FLOAT *rho, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
  void (*gga) (const struct XC(func_type) *p, int np, 
	       const FLOAT *rho, const FLOAT *sigma, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
	       FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3);
  void (*mgga)(const struct XC(func_type) *p, int np, 
	       const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	       FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2tau2, FLOAT *v2lapl2,
	       FLOAT *v2rhosigma, FLOAT *v2rhotau, FLOAT *v2rholapl, 
	       FLOAT *v2sigmatau, FLOAT *v2sigmalapl, FLOAT *v2taulapl);
} XC(func_info_type);

int XC(func_info_get_number)(const XC(func_info_type) *info);
int XC(func_info_get_kind)(const XC(func_info_type) *info);
char const *XC(func_info_get_name)(const XC(func_info_type) *info);
int XC(func_info_get_family)(const XC(func_info_type) *info);
int XC(func_info_get_flags)(const XC(func_info_type) *info);
const func_reference_type *XC(func_info_get_references)(const XC(func_info_type) *info, int number);
int XC(func_info_get_n_ext_params)(XC(func_info_type) *info);
char const *XC(func_info_get_ext_params_description)(XC(func_info_type) *info, int number);
double XC(func_info_get_ext_params_default_value)(XC(func_info_type) *info, int number);
  
struct XC(func_type){
  const XC(func_info_type) *info;       /* all the information concerning this functional */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int n_func_aux;                       /* how many auxiliary functions we need */
  struct XC(func_type) **func_aux;      /* most GGAs are based on a LDA or other GGAs  */
  FLOAT *mix_coef;                      /* coefficients for the mixing */

  FLOAT cam_omega;                      /* range-separation parameter for range-separated hybrids */
  FLOAT cam_alpha;                      /* fraction of Hartree-Fock exchange for normal or range-separated hybrids */
  FLOAT cam_beta;                       /* fraction of short-range exchange for range-separated hybrids */

  FLOAT nlc_b;                          /* Non-local correlation, b parameter */
  FLOAT nlc_C;                          /* Non-local correlation, C parameter */

  int n_rho, n_sigma, n_tau, n_lapl;    /* spin dimensions of the arrays */
  int n_zk;

  int n_vrho, n_vsigma, n_vtau, n_vlapl;

  int n_v2rho2, n_v2sigma2, n_v2tau2, n_v2lapl2,
    n_v2rhosigma, n_v2rhotau, n_v2rholapl, 
    n_v2sigmatau, n_v2sigmalapl, n_v2lapltau;

  int n_v3rho3, n_v3rho2sigma, n_v3rhosigma2, n_v3sigma3;

  void *params;                         /* this allows us to fix parameters in the functional */
};

typedef struct XC(func_type) XC(func_type);

/* functionals */
int   XC(functional_get_number)(const char *name);
char *XC(functional_get_name)(int number);
int   XC(family_from_id)(int id, int *family, int *number);
XC(func_type) *XC(func_alloc)();
int   XC(func_init)(XC(func_type) *p, int functional, int nspin);
void  XC(func_end)(XC(func_type) *p);
void  XC(func_free)(XC(func_type) *p);
const XC(func_info_type) *XC(func_get_info)(const XC(func_type) *p);
void  XC(func_set_ext_params)(XC(func_type) *p, double *ext_params);

#include "xc_funcs.h"
#include "xc_funcs_removed.h"
  
void XC(lda)        (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
void XC(lda_exc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk);
void XC(lda_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk, FLOAT *vrho);
void XC(lda_vxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *vrho);
void XC(lda_fxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *v2rho2);
void XC(lda_kxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *v3rho3);

void XC(gga)     (const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
		  FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3);
void XC(gga_exc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
		 FLOAT *zk);
void XC(gga_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma);
void XC(gga_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *vrho, FLOAT *vsigma);
void XC(gga_fxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
void XC(gga_kxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3);

void XC(gga_lb_modified)  (const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
			   FLOAT r, FLOAT *vrho);

FLOAT XC(gga_ak13_get_asymptotic) (FLOAT homo);

FLOAT XC(hyb_exx_coef)(const XC(func_type) *p);
void  XC(hyb_cam_coef)(const XC(func_type) *p, FLOAT *omega, FLOAT *alpha, FLOAT *beta);
void  XC(nlc_coef)(const XC(func_type) *p, FLOAT *nlc_b, FLOAT *nlc_C);

/* the meta-GGAs */
void XC(mgga)        (const XC(func_type) *p, int np,
		      const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
		      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau,
		      FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
		      FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
		      FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau);
void XC(mgga_exc)    (const XC(func_type) *p, int np,  
		      const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau, 
		      FLOAT *zk);
void XC(mgga_exc_vxc)(const XC(func_type) *p, int np, 
		      const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
		      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau);
void XC(mgga_vxc)    (const XC(func_type) *p, int np,
		      const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
		      FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau);
void XC(mgga_fxc)    (const XC(func_type) *p, int np, 
		      const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
		      FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
		      FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
		      FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau);

/* Functionals that are defined as mixtures of others */
void XC(mix_func)
  (const XC(func_type) *func, int np,
   const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl, const FLOAT *tau,
   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl, FLOAT *vtau,
   FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2lapl2, FLOAT *v2tau2,
   FLOAT *v2rhosigma, FLOAT *v2rholapl, FLOAT *v2rhotau, 
   FLOAT *v2sigmalapl, FLOAT *v2sigmatau, FLOAT *v2lapltau);
  
#include "xc_unconfig.h"

#ifdef __cplusplus
}
#endif

#endif
