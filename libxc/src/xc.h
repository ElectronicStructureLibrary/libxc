#ifndef _XC_H
#define _XC_H

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
  void (*lda) (const void *p, const double *rho, double *exc, double *vxc, double *fxc);
  void (*gga) (void *p, double *rho, double *sigma, double *exc, double *vrho, double *vsigma);
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
  double alpha;                     /* parameter for Xalpha functional */
} xc_lda_type;

int  xc_lda_init(xc_lda_type *p, int functional, int nspin);
void xc_lda_x_init(xc_lda_type *p, int nspin, int dim, int irel);
void xc_lda_c_xalpha_init(xc_lda_type *p, int nspin, int dim, double alpha);

void xc_lda(const xc_lda_type *p, const double *rho, double *exc, double *vxc, double *fxc, double *kxc);
void xc_lda_exc(const xc_lda_type *p, const double *rho, double *exc);
void xc_lda_vxc(const xc_lda_type *p, const double *rho, double *exc, double *vxc);
void xc_lda_fxc(const xc_lda_type *p, const double *rho, double *fxc);
void xc_lda_kxc(const xc_lda_type *p, const double *rho, double *kxc);

void xc_lda_sp(const xc_lda_type *p, const float *rho, float *exc, float *vxc, float *fxc, float *kxc);
void xc_lda_exc_sp(const xc_lda_type *p, const float *rho, float *exc);
void xc_lda_vxc_sp(const xc_lda_type *p, const float *rho, float *exc, float *vxc);
void xc_lda_fxc_sp(const xc_lda_type *p, const float *rho, float *fxc);
void xc_lda_kxc_sp(const xc_lda_type *p, const float *rho, float *kxc);


/* the GGAs */
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
void xc_gga_sp  (xc_gga_type *p, float *rho, float *grho,
		 float *e, float *dedd, float *dedgd);

void xc_gga_lb_init(xc_gga_type *p, int nspin, int modified, double threshold);
void xc_gga_lb     (xc_gga_type *p, double *rho, double *grho, double r, double ip, double qtot,
		    double *dedd);
void xc_gga_lb_sp  (xc_gga_type *p, float *rho, float *grho, float r, float ip, float qtot,
		    float *dedd);

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
void xc_mgga_sp  (xc_mgga_type *p, float *rho, float *grho, float *tau,
		  float *e, float *dedd, float *dedgd, float *dedtau);

/* the LCAs */

#define XC_LCA_OMC            301 /* Orestes, Marcasso & Capelle */
#define XC_LCA_LCH            302 /* Lee, Colwell & Handy        */


typedef struct{
  xc_func_info_type *info;  /* which functional did we chose   */
  int nspin;                /* XC_UNPOLARIZED or XC_POLARIZED  */

} xc_lca_type;

void xc_lca_init(xc_lca_type *p, int functional, int nspin);
void xc_lca     (xc_lca_type *p, double *rho, double *v, double *e, double *dedd, double *dedv);
void xc_lca_sp  (xc_lca_type *p, float *rho, float *v, float *e, float *dedd, float *dedv);

#ifdef __cplusplus
}
#endif

#endif

