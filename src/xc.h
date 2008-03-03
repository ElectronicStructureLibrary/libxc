/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef _XC_H
#define _XC_H

#include "xc_config.h"

#ifdef __cplusplus
extern "C" {
#endif
  
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
  void (*lda) (const void *p, const FLOAT *rho, FLOAT *exc, FLOAT *vxc, FLOAT *fxc);
  void (*gga) (void *p, FLOAT *rho, FLOAT *sigma, FLOAT *exc, FLOAT *vrho, FLOAT *vsigma);
} xc_func_info_type;


/* functionals */
int xc_family_from_id(int functional);
#include "xc_funcs.h"


/* the LDAs */
typedef struct struct_lda_type {
  const xc_func_info_type *info;    /* which functional did we chose   */
  int nspin;                        /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int relativistic;                 /* XC_RELATIVISTIC or XC_NON_RELATIVISTIC */
  int dim;
  
  struct struct_lda_type *lda_aux;  /* some LDAs are built on top of other LDAs */
  FLOAT alpha;                      /* parameter for Xalpha functional */
} xc_lda_type;

int  xc_lda_init(xc_lda_type *p, int functional, int nspin);
void xc_lda_x_init(xc_lda_type *p, int nspin, int dim, int irel);
void xc_lda_c_xalpha_init(xc_lda_type *p, int nspin, int dim, FLOAT alpha);

void xc_lda(const xc_lda_type *p, const FLOAT *rho, FLOAT *exc, FLOAT *vxc, FLOAT *fxc, FLOAT *kxc);
void xc_lda_exc(const xc_lda_type *p, const FLOAT *rho, FLOAT *exc);
void xc_lda_vxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *exc, FLOAT *vxc);
void xc_lda_fxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *fxc);
void xc_lda_kxc(const xc_lda_type *p, const FLOAT *rho, FLOAT *kxc);

/* the GGAs */
typedef struct xc_gga_type{
  const xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  xc_lda_type *lda_aux;           /* most GGAs are based on a LDA    */
  struct xc_gga_type **gga_aux;   /* and sometimes other GGAs */

  void *params;                   /* this allows to fix parameters in the functional */
} xc_gga_type;

int  xc_gga_init(xc_gga_type *p, int functional, int nspin);
void xc_gga_end (xc_gga_type *p);
void xc_gga     (xc_gga_type *p, FLOAT *rho, FLOAT *grho, FLOAT *e, FLOAT *dedd, FLOAT *dedgd);

void xc_gga_lb_set_params   (xc_gga_type *p, int modified, FLOAT threshold, FLOAT ip, FLOAT qtot);
void xc_gga_lb_modified     (xc_gga_type *p, FLOAT *rho, FLOAT *grho, FLOAT r, FLOAT *dedd);


/* the GGAs hybrids */
typedef struct xc_hyb_gga_type{
  const xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  xc_lda_type **lda_aux;          /* the LDA components of the hybrid */
  int           lda_n;            /* their number                     */
  FLOAT       *lda_coef;          /* and their coefficients           */

  xc_gga_type **gga_aux;          /* the GGA components               */
  int           gga_n;            /* their number                     */
  FLOAT       *gga_coef;          /* and their coefficients           */
  FLOAT        exx_coef;          /* the exact exchange coefficient   */

  void *params;                   /* this allows to fix parameters in the functional */
} xc_hyb_gga_type;

int  xc_hyb_gga_init(xc_hyb_gga_type *p, int functional, int nspin);
void xc_hyb_gga_end(xc_hyb_gga_type *p);
void xc_hyb_gga   (xc_hyb_gga_type *p, FLOAT *rho, FLOAT *sigma, FLOAT *e, FLOAT *vrho, FLOAT *vsigma);
FLOAT xc_hyb_gga_exx_coef   (xc_hyb_gga_type *p);


/* the meta-GGAs */

#define XC_MGGA_X_TPSS        201 /* Perdew, Tao, Staroverov & Scuseria exchange    */
#define XC_MGGA_C_TPSS        202 /* Perdew, Tao, Staroverov & Scuseria correlation */

typedef struct{
  const xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  xc_lda_type *lda_aux;     /* most meta-GGAs are based on a LDA    */
  xc_gga_type *gga_aux1;    /* or on a GGA                          */
  xc_gga_type *gga_aux2;    /* or on a GGA                          */

} xc_mgga_type;

void xc_mgga_init(xc_mgga_type *p, int functional, int nspin);
void xc_mgga_end (xc_mgga_type *p);
void xc_mgga     (xc_mgga_type *p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
		  FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau);

/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */

} xc_lca_type;

void xc_lca_init(xc_lca_type *p, int functional, int nspin);
void xc_lca     (xc_lca_type *p, FLOAT *rho, FLOAT *v, FLOAT *e, FLOAT *dedd, FLOAT *dedv);

#ifdef __cplusplus
}
#endif

#endif

