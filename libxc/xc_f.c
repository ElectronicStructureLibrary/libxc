#include <stdlib.h>
#include <stdio.h>

#include "xc.h"
#include "config.h"
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


void FC_FUNC_(xc_info_family, XC_INFO_FAMILY)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((func_type *)(*info))->family, s);
}

void FC_FUNC_(xc_info_ref, XC_INFO_REF)
     (void **info, int *n, STR_F_TYPE s STR_ARG1)
{
  func_type *func_p = (func_type *)(*info);
  if(func_p->refs[*n] == NULL)
    *n = -1;
  else{
    TO_F_STR1(func_p->refs[*n], s);
    (*n)++;
  }
}

/* LDAs */

void FC_FUNC_(xc_lda_init, XC_LDA_INIT)
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
  free(*p);
}

void FC_FUNC_(xc_lda, XC_LDA)
     (void **p, double *rho, double *e, double *v)
{
  lda((lda_type *)(*p), rho, e, v);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_lda_x_init, XC_LDA_X_INIT)
     (void **p, void **info, int *nspin, int *dim, int *rel)
{
  *p = malloc(sizeof(lda_type));
  lda_x_init((lda_type *)(*p), *nspin, *dim, *rel);
  *info = (void *)(((lda_type *) p)->func);
}

/* Slater's Xalpha */
void FC_FUNC_(xc_lda_c_xalpha_init, XC_LDA_C_XALPHA_INIT)
     (void **p, void **info, int *nspin, int *dim, int *rel, double *alpha)
{
  *p = malloc(sizeof(lda_type));
  lda_c_xalpha_init((lda_type *)(*p), *nspin, *dim, *rel, *alpha);
  *info = (void *)(((lda_type *) p)->func);
}


/* GGAs */

void FC_FUNC_(xc_gga_init, XC_GGA_INIT)
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
     (void **p, void **info,  int *nspin, int *modified, double *threshold)
{
  gga_type *gga_p;

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

