#include "config.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#ifndef _XC_H
#define _XC_H

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
  void (*lda) (void *p, double *rho, double *ec, double *vc, double *fc);
  void (*gga) (void *p, double *rho, double *sigma, double *ec, double *vc, double *vsigma);
} xc_func_info_type;


/* functionals */
int xc_family_from_id(int functional);


/* the LDAs */

#define XC_LDA_X                1   /* Exchange                     */
#define XC_LDA_C_WIGNER         2   /* Wigner parametrization       */
#define XC_LDA_C_RPA            3   /* Random Phase Approximation   */
#define XC_LDA_C_HL             4   /* Hedin & Lundqvist            */
#define XC_LDA_C_GL             5   /* Gunnarson & Lundqvist        */
#define XC_LDA_C_XALPHA         6   /* Slater's Xalpha              */
#define XC_LDA_C_VWN            7   /* Vosko, Wilk, & Nussair       */
#define XC_LDA_C_VWN_RPA        8   /* Vosko, Wilk, & Nussair (RPA) */
#define XC_LDA_C_PZ             9   /* Perdew & Zunger              */
#define XC_LDA_C_PZ_MOD        10   /* Perdew & Zunger (Modified)   */
#define XC_LDA_C_OB_PZ         11   /* Ortiz & Ballone (PZ)         */
#define XC_LDA_C_PW            12   /* Perdew & Wang                */
#define XC_LDA_C_OB_PW         13   /* Ortiz & Ballone (PW)         */
#define XC_LDA_C_AMGB          14   /* Attacalite et al             */

typedef struct struct_lda_type {
  const xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int relativistic;         /* XC_RELATIVISTIC or XC_NON_RELATIVISTIC */
  int dim;
  
  double alpha;             /* parameter for Xalpha functional */
} xc_lda_type;

int  xc_lda_init(xc_lda_type *p, int functional, int nspin);
void xc_lda_x_init(xc_lda_type *p, int nspin, int dim, int irel);
void xc_lda_c_xalpha_init(xc_lda_type *p, int nspin, int dim, double alpha);

void xc_lda(xc_lda_type *p, double *rho, double *ec, double *vc, double *fc);
void xc_lda_fxc(xc_lda_type *p, double *rho, double *fxc);
void xc_lda_kxc(xc_lda_type *p, double *rho, double *kxc);


/* the GGAs */

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange             */
#define XC_GGA_X_PBE_R        102 /* Perdew, Burke & Ernzerhof exchange (revised)   */
#define XC_GGA_X_B86          103 /* Becke 86 Xalfa,beta,gamma                      */
#define XC_GGA_X_B86_R        104 /* Becke 86 Xalfa,beta,gamma (reoptimized)        */
#define XC_GGA_X_B86_MGC      105 /* Becke 86 Xalfa,beta,gamma (with mod. grad. correction) */
#define XC_GGA_X_B88          106 /* Becke 88                                       */
#define XC_GGA_X_G96          107 /* Gill 96                                        */

#define XC_GGA_C_PBE          130 /* Perdew, Burke & Ernzerhof correlation          */
#define XC_GGA_C_LYP          131 /* Lee, Yang & Parr                               */

#define XC_GGA_XC_LB          160 /* van Leeuwen & Baerends                         */

typedef struct{
  const xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  xc_lda_type *lda_aux;     /* most GGAs are based on a LDA    */

  int modified;             /* parameters necessary to the lb functional */
  double threshold;
} xc_gga_type;

int  xc_gga_init(xc_gga_type *p, int functional, int nspin);
void xc_gga_end (xc_gga_type *p);
void xc_gga     (xc_gga_type *p, double *rho, double *grho,
		 double *e, double *dedd, double *dedgd);

void xc_gga_lb_init(xc_gga_type *p, int nspin, int modified, double threshold);
void xc_gga_lb     (xc_gga_type *p, double *rho, double *grho, double r, double ip, double qtot,
		 double *dedd);

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
void xc_mgga     (xc_mgga_type *p, double *rho, double *grho, double *tau,
		  double *e, double *dedd, double *dedgd, double *dedtau);

/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */

} xc_lca_type;

void xc_lca_init(xc_lca_type *p, int functional, int nspin);
void xc_lca     (xc_lca_type *p, double *rho, double *v, double *e, double *dedd, double *dedv);

#endif

