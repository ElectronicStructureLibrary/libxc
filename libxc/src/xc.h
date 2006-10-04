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

struct struct_lda_type;

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
} func_type;


/* functionals */
int family_from_id(int functional);


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

struct struct_lda_type {
  const func_type *func;      /* which functional did we chose   */
  int    nspin;         /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int    relativistic;  /* XC_RELATIVISTIC or XC_NON_RELATIVISTIC */
  int    dim;
  
  double alpha;         /* parameter for Xalpha functional */
};
typedef struct struct_lda_type lda_type;

int  lda_init(lda_type *p, int functional, int nspin);
void lda_x_init(lda_type *p, int nspin, int dim, int irel);
void lda_c_xalpha_init(lda_type *p, int nspin, int dim, double alpha);

void lda(lda_type *p, double *rho, double *ec, double *vc, double *fc);
void lda_fxc(lda_type *p, double *rho, double *fxc);
void lda_kxc(lda_type *p, double *rho, double *kxc);


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
  const func_type *func;       /* which functional did we chose   */
  int        nspin;      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  lda_type  *lda_aux;    /* most GGAs are based on a LDA    */

  int modified;          /* parameters necessary to the lb functional */
  double threshold;
} gga_type;

int  gga_init(gga_type *p, int functional, int nspin);
void gga_end (gga_type *p);
void gga     (gga_type *p, double *rho, double *grho,
	      double *e, double *dedd, double *dedgd);

void gga_lb_init(gga_type *p, int nspin, int modified, double threshold);
void gga_lb     (gga_type *p, double *rho, double *grho, double r, double ip, double qtot,
		 double *dedd);

/* the meta-GGAs */

#define XC_MGGA_X_TPSS        201 /* Perdew, Tao, Staroverov & Scuseria exchange    */
#define XC_MGGA_C_TPSS        202 /* Perdew, Tao, Staroverov & Scuseria correlation */

typedef struct{
  const func_type *func;       /* which functional did we chose   */
  int        nspin;      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  lda_type  *lda_aux;    /* most meta-GGAs are based on a LDA    */
  gga_type  *gga_aux1;   /* or on a GGA                          */
  gga_type  *gga_aux2;   /* or on a GGA                          */

} mgga_type;

void mgga_init(mgga_type *p, int functional, int nspin);
void mgga_end (mgga_type *p);
void mgga     (mgga_type *p, double *rho, double *grho, double *tau,
	       double *e, double *dedd, double *dedgd, double *dedtau);

/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  func_type *func;       /* which functional did we chose   */
  int        nspin;      /* XC_UNPOLARIZED or XC_POLARIZED  */

} lca_type;

void lca_init(lca_type *p, int functional, int nspin);
void lca     (lca_type *p, double *rho, double *v, double *e, double *dedd, double *dedv);

#endif

