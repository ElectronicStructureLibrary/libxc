#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "xc.h"
#include <string_f.h>

/* info */

int FC_FUNC_(xc_info_number, XC_INFO_NUMBER)
     (void **info)
{
  return ((func_type *)(*info))->number;
}


int FC_FUNC_(xc_info_kind, XC_INFO_KIND)
     (void **info)
{
  return ((func_type *)(*info))->kind;
}


void FC_FUNC_(xc_info_name, XC_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((func_type *)(*info))->name, s);
}


int FC_FUNC_(xc_info_family, XC_INFO_FAMILY)
     (void **info)
{
  return ((func_type *)(*info))->family;
}

void FC_FUNC_(xc_info_ref, XC_INFO_REF)
     (void **info, char **s, STR_F_TYPE ref_f STR_ARG1)
{
  char *c, ref[256]; /* hopefully no ref is longer than 256 characters ;) */
  func_type *func_p = (func_type *)(*info);

  if(*s == 0) *s = func_p->refs;

  if(*s == NULL || **s == '\0'){
    *s = (char *)(-1);
    return;
  }

  for(c=ref; **s!='\0' && **s!='\n'; (*s)++, c++)
    *c = **s;
  *c = '\0';
  if(**s=='\n') (*s)++;

  TO_F_STR1(ref, ref_f);
}

/* functionals */
int FC_FUNC_(xc_family_from_id, XC_FAMILY_FROM_ID)
  (int *functional)
{
  return family_from_id(*functional);
}


/* LDAs */

void FC_FUNC_(xc_lda_init_, XC_LDA_INIT_)
     (void **p, void **info, int *functional, int *nspin)
{
  lda_type *lda_p;
  
  *p = malloc(sizeof(lda_type));
  lda_p = (lda_type *)(*p);
  lda_init(lda_p, *functional, *nspin);
  *info = (void *)(lda_p->func);
}

void FC_FUNC_(xc_lda_end, XC_LDA_END)
     (void **p)
{
}

void FC_FUNC_(xc_lda, XC_LDA)
     (void **p, double *rho, double *e, double *v)
{
  lda((lda_type *)(*p), rho, e, v, NULL);
}

void FC_FUNC_(xc_lda_fxc, XC_LDA_FXC)
     (void **p, double *rho, double *fxc)
{
  lda_fxc((lda_type *)(*p), rho, fxc);
}

void FC_FUNC_(xc_lda_kxc, XC_LDA_KXC)
     (void **p, double *rho, double *kxc)
{
  lda_kxc((lda_type *)(*p), rho, kxc);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_lda_x_init, XC_LDA_X_INIT)
     (void **p, void **info, int *functional, int *nspin, int *dim, int *irel)
{
  lda_type *lda_p;

  assert(*functional == XC_LDA_X);

  *p = malloc(sizeof(lda_type));
  lda_p = (lda_type *)(*p);
  lda_x_init(lda_p, *nspin, *dim, *irel);
  *info = (void *)(lda_p->func);
}

/* Slater's Xalpha */
void FC_FUNC_(xc_lda_c_xalpha_init, XC_LDA_C_XALPHA_INIT)
     (void **p, void **info, int *functional, int *nspin, int *dim, double *alpha)
{
  lda_type *lda_p;

  assert(*functional == XC_LDA_C_XALPHA);

  *p = malloc(sizeof(lda_type));
  lda_p = (lda_type *)(*p);
  lda_c_xalpha_init(lda_p, *nspin, *dim, *alpha);
  *info = (void *)(lda_p->func);
}


/* GGAs */

void FC_FUNC_(xc_gga_init_, XC_GGA_INIT_)
     (void **p, void **info, int *functional, int *nspin)
{
  gga_type *gga_p;

  *p = malloc(sizeof(gga_type));
  gga_p = (gga_type *)(*p);
  gga_init(gga_p, *functional, *nspin);
  *info = (void *)(gga_p->func);
}

void FC_FUNC_(xc_gga_end, XC_GGA_END)
     (void **p)
{
  gga_end((gga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_gga, XC_GGA)
     (void **p, double *rho, double *grho, 
      double *e, double *dedd, double *dedgd)
{
  gga((gga_type *)(*p), rho, grho, e, dedd, dedgd);
}


/* the van Leeuwen & Baerends functional is special */
void FC_FUNC_(xc_gga_lb_init, XC_GGA_LB_INIT)
     (void **p, void **info,  int *functional, int *nspin, int *modified, double *threshold)
{
  gga_type *gga_p;

  assert(*functional == XC_GGA_XC_LB);

  *p = malloc(sizeof(gga_type));
  gga_p = (gga_type *)(*p);
  gga_lb_init(gga_p, *nspin, *modified, *threshold);
  *info = (void *)(gga_p->func);
}

void FC_FUNC_(xc_gga_lb, XC_GGA_LB)
     (void **p, double *rho, double *grho, double *r, double *ip, double *qtot,
      double *dedd)
{
  gga_lb((gga_type *)(*p), rho, grho, *r, *ip, *qtot, dedd);
}


/* meta-GGAs */

void FC_FUNC_(xc_mgga_init, XC_MGGA_INIT)
     (void **p, void **info, int *functional, int *nspin)
{
  mgga_type *mgga_p;

  *p = malloc(sizeof(mgga_type));
  mgga_p = (mgga_type *)(*p);
  mgga_init(mgga_p, *functional, *nspin);
  *info = (void *)(mgga_p->func);
}

void FC_FUNC_(xc_mgga_end, XC_MGGA_END)
     (void **p)
{
  mgga_end((mgga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_mgga, XC_MGGA)
  (void **p, double *rho, double *grho, double *tau,
   double *e, double *dedd, double *dedgd, double *dedtau)
{
  mgga((mgga_type *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
}


/* LCAs */

void FC_FUNC_(xc_lca_init, XC_LCA_INIT)
     (void **p, void **info, int *functional, int *nspin)
{
  lca_type *lca_p;

  *p = malloc(sizeof(lca_type));
  lca_p = (lca_type *)(*p);
  lca_init(lca_p, *functional, *nspin);
  *info = (void *)(lca_p->func);
}

void FC_FUNC_(xc_lca_end, XC_LCA_END)
     (void **p)
{
  free(*p);
}

void FC_FUNC_(xc_lca, XC_LCA)
     (void **p, double *rho, double *v, 
      double *e, double *dedd, double *dedv)
{
  lca((lca_type *)(*p), rho, v, e, dedd, dedv);
}
