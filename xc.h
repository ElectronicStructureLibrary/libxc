#ifndef _XC_H
#define _XC_H


#define XC_UNPOLARIZED          1
#define XC_POLARIZED            2

#define XC_NON_RELATIVISTIC     0
#define XC_RELATIVISTIC         1

#define XC_EXCHANGE             0
#define XC_CORRELATION          1
#define XC_EXCHANGE_CORRELATION 2

typedef struct{
  int   number; /* indentifier number */
  int   kind;   /* XC_EXCHANGE or XC_CORRELATION */

  char *name;   /* name of the functional, e.g. PBE */
  char *family; /* type of the functional, e.g. GGA */
  char *refs;  /* references                       */
}func_type;

/* the LDAs */

#define XC_LDA_X                1   /* Exchange                   */
#define XC_LDA_C_WIGNER         2   /* Wigner parametrization     */
#define XC_LDA_C_RPA            3   /* Random Phase Approximation */
#define XC_LDA_C_HL             4   /* Hedin & Lundqvist          */
#define XC_LDA_C_GL             5   /* Gunnarson & Lundqvist      */
#define XC_LDA_C_XALPHA         6   /* Slater's Xalpha            */
#define XC_LDA_C_VWN            7   /* Vosko, Wilk, & Nussair     */
#define XC_LDA_C_PZ             8   /* Perdew & Zunger            */
#define XC_LDA_C_OB_PZ          9   /* Ortiz & Ballone (PZ)       */
#define XC_LDA_C_PW            10   /* Perdew & Wang              */
#define XC_LDA_C_OB_PW         11   /* Ortiz & Ballone (PW)       */
#define XC_LDA_C_LYP           12   /* Lee, Yang, & Parr LDA      */
#define XC_LDA_C_AMGB          13   /* Attacalite et al           */

typedef struct{
  func_type *func;      /* which functional did we chose   */
  int    nspin;         /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  int    relativistic;  /* not used for the moment         */
  int    dim;
  
  double alpha;         /* parameter for Xalpha functional */
} lda_type;

void lda_init(lda_type *p, int functional, int nspin);
void lda_x_init(lda_type *p, int nspin, int dim);
void lda_c_xalpha_init(lda_type *p, int nspin, int dim, double alpha);

void lda(lda_type *p, double *rho, double *ec, double *vc);
void lda_fxc(lda_type *p, double *rho, double *fxc);


/* the GGAs */

#define XC_GGA_X_PBE          101 /* Perdew, Burke & Ernzerhof exchange    */
#define XC_GGA_C_PBE          102 /* Perdew, Burke & Ernzerhof correlation */
#define XC_GGA_XC_LB          103 /* van Leeuwen & Baerends                */

typedef struct{
  func_type *func;       /* which functional did we chose   */
  int        nspin;      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  lda_type  *lda_aux;    /* most GGAs are based on a LDA    */

  int modified;          /* parameters necessary to the lb functional */
  double threshold;
} gga_type;

void gga_init(gga_type *p, int functional, int nspin);
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
  func_type *func;       /* which functional did we chose   */
  int        nspin;      /* XC_UNPOLARIZED or XC_POLARIZED  */
  
  lda_type  *lda_aux;    /* most meta-GGAs are based on a LDA    */
  gga_type  *gga_aux1;   /* or on a GGA                          */
  gga_type  *gga_aux2;   /* or on a GGA                          */

} mgga_type;

void mgga_init(mgga_type *p, int functional, int nspin);
void mgga_end (mgga_type *p);
void mgga     (mgga_type *p, double *rho, double *grho, double *tau,
	       double *e, double *dedd, double *dedgd, double *dedtau);

#endif
