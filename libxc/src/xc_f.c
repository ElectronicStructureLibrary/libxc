#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "xc.h"
#include <string_f.h>

/* info */

int FC_FUNC_(xc_f90_info_number, XC_F90_INFO_NUMBER)
     (void **info)
{
  return ((xc_func_info_type *)(*info))->number;
}


int FC_FUNC_(xc_f90_info_kind, XC_F90_INFO_KIND)
     (void **info)
{
  return ((xc_func_info_type *)(*info))->kind;
}


void FC_FUNC_(xc_f90_info_name, XC_F90_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((xc_func_info_type *)(*info))->name, s);
}


int FC_FUNC_(xc_f90_info_family, XC_F90_INFO_FAMILY)
     (void **info)
{
  return ((xc_func_info_type *)(*info))->family;
}

void FC_FUNC_(xc_f90_info_ref, XC_F90_INFO_REF)
     (void **info, char **s, STR_F_TYPE ref_f STR_ARG1)
{
  char *c, ref[256]; /* hopefully no ref is longer than 256 characters ;) */
  xc_func_info_type *func_p = (xc_func_info_type *)(*info);

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
int FC_FUNC_(xc_f90_family_from_id, XC_F90_FAMILY_FROM_ID)
  (int *functional)
{
  return xc_family_from_id(*functional);
}


/* LDAs */

void FC_FUNC_(xc_f90_lda_init_, XC_F90_LDA_INIT_)
     (void **p, void **info, int *functional, int *nspin)
{
  xc_lda_type *lda_p;
  
  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_init(lda_p, *functional, *nspin);
  *info = (void *)(lda_p->info);
}

void FC_FUNC_(xc_f90_lda_end, XC_F90_LDA_END)
     (void **p)
{
}

void FC_FUNC_(xc_f90_lda, XC_F90_LDA)
     (void **p, double *rho, double *e, double *v)
{
  xc_lda((xc_lda_type *)(*p), rho, e, v, NULL);
}

void FC_FUNC_(xc_f90_lda_fxc, XC_F90_LDA_FXC)
     (void **p, double *rho, double *fxc)
{
  xc_lda_fxc((xc_lda_type *)(*p), rho, fxc);
}

void FC_FUNC_(xc_f90_lda_kxc, XC_F90_LDA_KXC)
     (void **p, double *rho, double *kxc)
{
  xc_lda_kxc((xc_lda_type *)(*p), rho, kxc);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_f90_lda_x_init, XC_F90_LDA_X_INIT)
     (void **p, void **info, int *functional, int *nspin, int *dim, int *irel)
{
  xc_lda_type *lda_p;

  assert(*functional == XC_LDA_X);

  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_x_init(lda_p, *nspin, *dim, *irel);
  *info = (void *)(lda_p->info);
}

/* Slater's Xalpha */
void FC_FUNC_(xc_f90_lda_c_xalpha_init, XC_F90_LDA_C_XALPHA_INIT)
     (void **p, void **info, int *functional, int *nspin, int *dim, double *alpha)
{
  xc_lda_type *lda_p;

  assert(*functional == XC_LDA_C_XALPHA);

  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_c_xalpha_init(lda_p, *nspin, *dim, *alpha);
  *info = (void *)(lda_p->info);
}


/* GGAs */

void FC_FUNC_(xc_f90_gga_init_, XC_F90_GGA_INIT_)
     (void **p, void **info, int *functional, int *nspin)
{
  xc_gga_type *gga_p;

  *p = malloc(sizeof(xc_gga_type));
  gga_p = (xc_gga_type *)(*p);
  xc_gga_init(gga_p, *functional, *nspin);
  *info = (void *)(gga_p->info);
}

void FC_FUNC_(xc_f90_gga_end, XC_F90_GGA_END)
     (void **p)
{
  xc_gga_end((xc_gga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_f90_gga, XC_F90_GGA)
     (void **p, double *rho, double *grho, 
      double *e, double *dedd, double *dedgd)
{
  xc_gga((xc_gga_type *)(*p), rho, grho, e, dedd, dedgd);
}


/* the van Leeuwen & Baerends functional is special */
void FC_FUNC_(xc_f90_gga_lb_init, XC_F90_GGA_LB_INIT)
     (void **p, void **info,  int *functional, int *nspin, int *modified, double *threshold)
{
  xc_gga_type *gga_p;

  assert(*functional == XC_GGA_XC_LB);

  *p = malloc(sizeof(xc_gga_type));
  gga_p = (xc_gga_type *)(*p);
  xc_gga_lb_init(gga_p, *nspin, *modified, *threshold);
  *info = (void *)(gga_p->info);
}

void FC_FUNC_(xc_f90_gga_lb, XC_F90_GGA_LB)
     (void **p, double *rho, double *grho, double *r, double *ip, double *qtot,
      double *dedd)
{
  xc_gga_lb((xc_gga_type *)(*p), rho, grho, *r, *ip, *qtot, dedd);
}


/* meta-GGAs */

void FC_FUNC_(xc_f90_mgga_init, XC_F90_MGGA_INIT)
     (void **p, void **info, int *functional, int *nspin)
{
  xc_mgga_type *mgga_p;

  *p = malloc(sizeof(xc_mgga_type));
  mgga_p = (xc_mgga_type *)(*p);
  xc_mgga_init(mgga_p, *functional, *nspin);
  *info = (void *)(mgga_p->info);
}

void FC_FUNC_(xc_f90_mgga_end, XC_F90_MGGA_END)
     (void **p)
{
  xc_mgga_end((xc_mgga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_f90_mgga, XC_F90_MGGA)
  (void **p, double *rho, double *grho, double *tau,
   double *e, double *dedd, double *dedgd, double *dedtau)
{
  xc_mgga((xc_mgga_type *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
}


/* LCAs */

void FC_FUNC_(xc_f90_lca_init, XC_F90_LCA_INIT)
     (void **p, void **info, int *functional, int *nspin)
{
  xc_lca_type *lca_p;

  *p = malloc(sizeof(xc_lca_type));
  lca_p = (xc_lca_type *)(*p);
  xc_lca_init(lca_p, *functional, *nspin);
  *info = (void *)(lca_p->info);
}

void FC_FUNC_(xc_f90_lca_end, XC_F90_LCA_END)
     (void **p)
{
  free(*p);
}

void FC_FUNC_(xc_f90_lca, XC_F90_LCA)
     (void **p, double *rho, double *v, 
      double *e, double *dedd, double *dedv)
{
  xc_lca((xc_lca_type *)(*p), rho, v, e, dedd, dedv);
}
