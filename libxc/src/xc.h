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

#include "xc_config.h"
  
#define XC_UNPOLARIZED          1
#define XC_POLARIZED            2

#define XC_NON_RELATIVISTIC     0
#define XC_RELATIVISTIC         1

#define XC_EXCHANGE             0
#define XC_CORRELATION          1
#define XC_EXCHANGE_CORRELATION 2

#define XC_FAMILY_UNKNOWN      -1
#define XC_FAMILY_LDA           1
#define XC_FAMILY_GGA           2
#define XC_FAMILY_MGGA          4
#define XC_FAMILY_LCA           8
#define XC_FAMILY_OEP          16
#define XC_FAMILY_HYB_GGA      32

#define XC_PROVIDES_EXC         1
#define XC_PROVIDES_VXC         2
#define XC_PROVIDES_FXC         4
#define XC_PROVIDES_KXC         8

typedef struct{
  int   number;   /* indentifier number */
  int   kind;     /* XC_EXCHANGE or XC_CORRELATION */

  char *name;     /* name of the functional, e.g. "PBE" */
  int   family;   /* type of the functional, e.g. XC_FAMILY_GGA */
  char *refs;     /* references                       */

  int   provides; /* what the functional provides, e.g. XC_PROVIDES_EXC | XC_PROVIDES_VXC */

  void (*init)(void *p);
  void (*end) (void *p);
  void (*lda) (const void *p, int np, const FLOAT *rho, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
  void (*gga) (const void *p, int np, const FLOAT *rho, const FLOAT *sigma, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
  void (*mgga)(const void *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);
} XC(func_info_type);


struct XC(struct_lda_type);
struct XC(struct_gga_type);
struct XC(struct_mgga_type);

typedef struct XC(struct_func_type){
  const XC(func_info_type) *info;       /* all the information concerning this functional */
  int nspin;                            /* this is a copy from the underlying functional */

  struct XC(struct_lda_type)  *lda;
  struct XC(struct_gga_type)  *gga;
  struct XC(struct_mgga_type) *mgga;
} XC(func_type);


/* functionals */
int  XC(family_from_id)(int id, int *family, int *number);
int  XC(func_init)(XC(func_type) *p, int functional, int nspin);
void XC(func_end)(XC(func_type) *p);

#include "xc_funcs.h"

/* the LDAs */
typedef struct XC(struct_lda_type) {
  const XC(func_info_type) *info;       /* all the information concerning this functional */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */

  int func;                             /* Shortcut in case of several functionals sharing the same interface */
  int n_rho, n_zk, n_vrho, n_v2rho2, n_v3rho3; /* spin dimensions of arguments */

  void *params;                         /* this allows us to fix parameters in the functional */
} XC(lda_type);

int  XC(lda_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin);
void XC(lda_end) (XC(func_type) *p);

void XC(lda)        (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
void XC(lda_exc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk);
void XC(lda_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk, FLOAT *vrho);
void XC(lda_vxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *vrho);
void XC(lda_fxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *v2rho2);
void XC(lda_kxc)    (const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *v3rho3);

void XC(lda_x_1d_set_params)     (XC(func_type) *p, int interaction, FLOAT bb);
void XC(lda_c_1d_csc_set_params) (XC(func_type) *p, int interaction, FLOAT bb);
void XC(lda_c_xalpha_set_params) (XC(func_type) *p, FLOAT alpha);
void XC(lda_x_set_params)        (XC(func_type) *p, int relativistic);
void XC(lda_c_2d_prm_set_params) (XC(func_type) *p, FLOAT N);
void XC(lda_c_vwn_set_params)    (XC(func_type) *p, int spin_interpolation);


/* the GGAs */
typedef struct XC(struct_gga_type){
  const XC(func_info_type) *info;       /* which functional did we choose   */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED   */
  
  int n_func_aux;                       /* how many auxiliary functions we need */
  XC(func_type) **func_aux;             /* most GGAs are based on a LDA or other GGAs  */
  FLOAT *mix_coef;                      /* coefficients for the mixing */

  FLOAT exx_coef;                       /* the Hartree-Fock mixing parameter for the hybrids */

  int func;                             /* Shortcut in case of several functionals sharing the same interface */
  int n_rho, n_zk, n_vrho, n_v2rho2;    /* spin dimensions of arguments */
  int n_sigma, n_vsigma, n_v2rhosigma, n_v2sigma2;

  void *params;                         /* this allows us to fix parameters in the functional */
} XC(gga_type);

int  XC(gga_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin);
void XC(gga_end) (XC(func_type) *p);
void XC(gga)     (const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
void XC(gga_exc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
		 FLOAT *zk);
void XC(gga_exc_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma);
void XC(gga_vxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *vrho, FLOAT *vsigma);
void XC(gga_fxc)(const XC(func_type) *p, int np, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);

void XC(gga_lb_modified)  (const XC(gga_type) *p, int np, const FLOAT *rho, const FLOAT *sigma, 
			   FLOAT r, FLOAT *vrho);

void XC(gga_x_b88_set_params)(XC(func_type) *p, FLOAT beta, FLOAT gamma);
void XC(gga_c_lyp_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d);
void XC(gga_lb_set_params)   (XC(func_type) *p, int modified, FLOAT threshold, FLOAT ip, FLOAT qtot);

FLOAT XC(hyb_gga_exx_coef)(XC(gga_type) *p);


/* the meta-GGAs */
typedef struct XC(struct_mgga_type){
  const XC(func_info_type) *info;       /* which functional did we choose   */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int n_func_aux;                       /* how many auxiliary functions we need */
  XC(func_type) **func_aux;             /* most GGAs are based on a LDA or other GGAs  */

  int func;                             /* Shortcut in case of several functionals sharing the same interface */

  void *params;                         /* this allows us to fix parameters in the functional */
} XC(mgga_type);

int  XC(mgga_init)(XC(func_type) *p, const XC(func_info_type) *info, int nspin);
void XC(mgga_end) (XC(func_type) *p);
void XC(mgga)(const XC(func_type) *p, const FLOAT *rho, 
	      const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);
void XC(mgga_exc)(const XC(func_type) *p, const FLOAT *rho, 
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau, 
		  FLOAT *zk);
void XC(mgga_exc_vxc)(const XC(func_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau);
void XC(mgga_vxc)(const XC(func_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau);
void XC(mgga_fxc)(const XC(func_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);

void XC(mgga_x_tb09_set_params)(XC(func_type) *p, FLOAT c);

/* Functionals that are defined as mixtures of others */
void XC(mix_func)(const XC(func_type) *dest_func, int n_func_aux, XC(func_type) **func_aux, FLOAT *mix_coef,
		  int np, const FLOAT *rho, const FLOAT *sigma,
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
  
/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  XC(func_info_type) *info;  /* which functional did we choose   */
  int nspin;                 /* XC_UNPOLARIZED or XC_POLARIZED  */

} XC(lca_type);

void XC(lca_init)(XC(lca_type) *p, int functional, int nspin);
void XC(lca)     (XC(lca_type) *p, FLOAT *rho, FLOAT *v, FLOAT *e, FLOAT *dedd, FLOAT *dedv);

#ifdef __cplusplus
}
#endif

#endif

