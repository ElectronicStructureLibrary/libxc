#ifndef _XC_H
#define _XC_H


#define XC_UNPOLARIZED          1
#define XC_POLARIZED            2

#define XC_NON_RELATIVISTIC     0
#define XC_RELATIVISTIC         1

/* the LDA */

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
	int    functional;    /* which functional did we chose   */
	int    nspin;         /* XC_UNPOLARIZED or XC_POLARIZED  */

  int    relativistic;  /* necessary for the exchange      */
	int    dim;
	
	double alpha;         /* parameter for Xalpha functional */
} lda_type;

void lda_init(lda_type *p, int nspin, int functional);
void lda_x_init(lda_type *p, int nspin, int dim, int rel);
void lda_c_xalpha_init(lda_type *p, int nspin, int dim, int rel, double alpha);

void lda(lda_type *p, double *rho, double *ec, double *vc);


#endif
