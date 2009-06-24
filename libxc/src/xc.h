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
  void (*lda) (const void *p, const FLOAT *rho, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
  void (*gga) (const void *p, const FLOAT *rho, const FLOAT *sigma, 
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
  void (*mgga)(const void *p, const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);
} XC(func_info_type);


/* functionals */
int XC(family_from_id)(int functional);
#include "xc_funcs.h"

struct XC(struct_mix_func_type);

/* the LDAs */
typedef struct XC(struct_lda_type) {
  const XC(func_info_type) *info;       /* which functional did we chose   */
  int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int relativistic;                     /* XC_RELATIVISTIC or XC_NON_RELATIVISTIC */
  int dim;
  
  struct XC(struct_lda_type) *lda_aux;  /* some LDAs are built on top of other LDAs */
  struct XC(struct_mix_func_type) *mix; /* others are mixtures of LDAs */

  void *params;                         /* this allows to fix parameters in the functional */
  FLOAT alpha;                          /* parameter for Xalpha functional (to disappear) */
} XC(lda_type);

int  XC(lda_init)(XC(lda_type) *p, int functional, int nspin);
void XC(lda_end) (XC(lda_type) *p);
void XC(lda_x_init)(XC(lda_type) *p, int nspin, int dim, int irel);
void XC(lda_c_xalpha_init)(XC(lda_type) *p, int nspin, int dim, FLOAT alpha);

void XC(lda)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3);
void XC(lda_exc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *zk);
void XC(lda_exc_vxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *zk, FLOAT *vrho);
void XC(lda_vxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *vrho);
void XC(lda_fxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *v2rho2);
void XC(lda_kxc)(const XC(lda_type) *p, const FLOAT *rho, FLOAT *v3rho3);

void XC(lda_c_1d_csc_set_params)  (XC(lda_type) *p, FLOAT bb);
void XC(lda_c_2d_prm_set_params)(XC(lda_type) *p, FLOAT N);
void XC(lda_c_vwn_set_params)     (const XC(lda_type) *p, int spin_interpolation);


/* the GGAs */
typedef struct XC(struct_gga_type){
  const XC(func_info_type) *info;         /* which functional did we chose   */
  int nspin;                              /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  XC(lda_type) *lda_aux;                  /* most GGAs are based on a LDA    */

  struct XC(struct_gga_type) **gga_aux;   /* and sometimes other GGAs */
  struct XC(struct_mix_func_type) *mix;   /* others are mixtures of LDAs and GGAs */

  void *params;                           /* this allows to fix parameters in the functional */
} XC(gga_type);

int  XC(gga_init)(XC(gga_type) *p, int functional, int nspin);
void XC(gga_end) (XC(gga_type) *p);
void XC(gga)     (const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
void XC(gga_exc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		 FLOAT *zk);
void XC(gga_exc_vxc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma);
void XC(gga_vxc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *vrho, FLOAT *vsigma);
void XC(gga_fxc)(const XC(gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);

void XC(gga_lb_set_params)   (XC(gga_type) *p, int modified, FLOAT threshold, FLOAT ip, FLOAT qtot);
void XC(gga_lb_modified)     (XC(gga_type) *p, FLOAT *rho, FLOAT *grho, FLOAT r, FLOAT *vrho);


/* the GGAs hybrids */
typedef struct XC(struct_hyb_gga_type){
  const XC(func_info_type) *info;  /* which functional did we chose   */
  int nspin;                       /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  struct XC(struct_mix_func_type) *mix; /* these functionals are usually 
					   mixtures of LDAs and GGAs */

  FLOAT exx_coef;                  /* the Hartree-Fock mixing parameter */
  void *params;                    /* this allows to fix parameters in the functional */
} XC(hyb_gga_type);

int  XC(hyb_gga_init)(XC(hyb_gga_type) *p, int functional, int nspin);
void XC(hyb_gga_end)(XC(hyb_gga_type) *p);
void XC(hyb_gga)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
		 FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
void XC(hyb_gga_exc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma, 
		     FLOAT *zk);
void XC(hyb_gga_exc_vxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
			 FLOAT *zk, FLOAT *vrho, FLOAT *vsigma);
void XC(hyb_gga_vxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		     FLOAT *vrho, FLOAT *vsigma);
void XC(hyb_gga_fxc)(const XC(hyb_gga_type) *p, const FLOAT *rho, const FLOAT *sigma,
		     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);
FLOAT XC(hyb_gga_exx_coef)(XC(hyb_gga_type) *p);


/* the meta-GGAs */
typedef struct{
  const XC(func_info_type) *info;  /* which functional did we chose   */
  int nspin;                       /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  XC(lda_type) *lda_aux;           /* most meta-GGAs are based on a LDA    */
  XC(gga_type) *gga_aux1;          /* or on a GGA                          */
  XC(gga_type) *gga_aux2;          /* or on a GGA                          */

  void *params;                    /* this allows to fix parameters in the functional */
} XC(mgga_type);

int  XC(mgga_init)(XC(mgga_type) *p, int functional, int nspin);
void XC(mgga_end) (XC(mgga_type) *p);
void XC(mgga)(const XC(mgga_type) *p, const FLOAT *rho, 
	      const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
	      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
	      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);
void XC(mgga_exc)(const XC(mgga_type) *p, const FLOAT *rho, 
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau, 
		  FLOAT *zk);
void XC(mgga_exc_vxc)(const XC(mgga_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau);
void XC(mgga_vxc)(const XC(mgga_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau);
void XC(mgga_fxc)(const XC(mgga_type) *p, const FLOAT *rho,
		  const FLOAT *sigma, const FLOAT *lapl_rho, const FLOAT *tau,
		  FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2);

/* Functionals that are defined as mixtures of others */
typedef struct XC(struct_mix_func_type){
  int level; /* LDA, GGA, MGGA */
  int nspin; /* POLARIZED, UNPOLARIZED */

  int            lda_n;
  XC(lda_type)  *lda_mix;
  FLOAT         *lda_coef;

  int            gga_n;
  XC(gga_type)  *gga_mix;
  FLOAT         *gga_coef;

  int            mgga_n;
  XC(mgga_type) *mgga_mix;
  FLOAT         *mgga_coef;
} XC(mix_func_type);

void XC(mix_func_init)(XC(mix_func_type) *p, int level, int nspin);
void XC(mix_func_alloc)(XC(mix_func_type) *p);
void XC(mix_func_free)(XC(mix_func_type) *p);
void XC(mix_func)(XC(mix_func_type) *p, const FLOAT *rho, const FLOAT *sigma,
		  FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2);

/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  XC(func_info_type) *info;  /* which functional did we chose   */
  int nspin;                 /* XC_UNPOLARIZED or XC_POLARIZED  */

} XC(lca_type);

void XC(lca_init)(XC(lca_type) *p, int functional, int nspin);
void XC(lca)     (XC(lca_type) *p, FLOAT *rho, FLOAT *v, FLOAT *e, FLOAT *dedd, FLOAT *dedv);

#ifdef __cplusplus
}
#endif

#endif

