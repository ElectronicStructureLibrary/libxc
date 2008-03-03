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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "config.h"

#ifdef HAVE_FORTRAN

#include "xc.h"
#include <string_f.h>

/* info */

CC_FORTRAN_INT FC_FUNC_(xc_f90_info_number, XC_F90_INFO_NUMBER)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->number;
}


CC_FORTRAN_INT FC_FUNC_(xc_f90_info_kind, XC_F90_INFO_KIND)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->kind;
}


void FC_FUNC_(xc_f90_info_name, XC_F90_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((xc_func_info_type *)(*info))->name, s);
}


CC_FORTRAN_INT  FC_FUNC_(xc_f90_info_family, XC_F90_INFO_FAMILY)
     (void **info)
{
  return (CC_FORTRAN_INT) ((xc_func_info_type *)(*info))->family;
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
CC_FORTRAN_INT  FC_FUNC_(xc_f90_family_from_id, XC_F90_FAMILY_FROM_ID)
  (CC_FORTRAN_INT  *functional)
{
  return (CC_FORTRAN_INT) xc_family_from_id((int) (*functional));
}


/* LDAs */

/* Standard initialization */
void FC_FUNC_(xc_f90_lda_init_, XC_F90_LDA_INIT_)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_lda_type *lda_p;
  
  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_init(lda_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(lda_p->info);
}

void FC_FUNC_(xc_f90_lda_end, XC_F90_LDA_END)
     (void **p)
{
  /* xc_lda_end does not exist */
  free(*p);
  *p = NULL;
}

void FC_FUNC_(xc_f90_lda, XC_F90_LDA)
  (void **p, FLOAT *rho, FLOAT *exc, FLOAT *vxc, FLOAT *fxc, FLOAT *kxc)
{
  xc_lda((xc_lda_type *)(*p), rho, exc, vxc, fxc, kxc);
}

void FC_FUNC_(xc_f90_lda_vxc, XC_F90_LDA_VXC)
     (void **p, FLOAT *rho, FLOAT *e, FLOAT *v)
{
  xc_lda_vxc((xc_lda_type *)(*p), rho, e, v);
}

void FC_FUNC_(xc_f90_lda_fxc, XC_F90_LDA_FXC)
     (void **p, FLOAT *rho, FLOAT *fxc)
{
  xc_lda_fxc((xc_lda_type *)(*p), rho, fxc);
}

void FC_FUNC_(xc_f90_lda_kxc, XC_F90_LDA_KXC)
     (void **p, FLOAT *rho, FLOAT *kxc)
{
  xc_lda_kxc((xc_lda_type *)(*p), rho, kxc);
}


/* Now come some special initializations */

/* exchange in the LDA */
void FC_FUNC_(xc_f90_lda_x_init, XC_F90_LDA_X_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, CC_FORTRAN_INT *irel)
{
  xc_lda_type *lda_p;

  assert(*functional == XC_LDA_X);

  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_x_init(lda_p, (int) (*nspin), (int) (*dim), (int) (*irel));
  *info = (void *)(lda_p->info);
}

/* Slater's Xalpha */
void FC_FUNC_(xc_f90_lda_c_xalpha_init, XC_F90_LDA_C_XALPHA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, FLOAT *alpha)
{
  xc_lda_type *lda_p;

  assert((int) (*functional) == XC_LDA_C_XALPHA);

  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_c_xalpha_init(lda_p, (int) (*nspin), (int) (*dim), *alpha);
  *info = (void *)(lda_p->info);
}


/* GGAs */

void FC_FUNC_(xc_f90_gga_init, XC_F90_GGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_gga_type *gga_p;

  *p = malloc(sizeof(xc_gga_type));
  gga_p = (xc_gga_type *)(*p);
  xc_gga_init(gga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(gga_p->info);
}

void FC_FUNC_(xc_f90_gga_end, XC_F90_GGA_END)
     (void **p)
{
  xc_gga_end((xc_gga_type *)(*p));
  free(*p);
  *p = NULL;
}

void FC_FUNC_(xc_f90_gga, XC_F90_GGA)
     (void **p, FLOAT *rho, FLOAT *grho, 
      FLOAT *e, FLOAT *dedd, FLOAT *dedgd)
{
  xc_gga((xc_gga_type *)(*p), rho, grho, e, dedd, dedgd);
}


/* the van Leeuwen & Baerends functional is special */
void FC_FUNC_(xc_f90_gga_lb_set_params, XC_F90_GGA_LB_SET_PARAMS)
  (void **p, CC_FORTRAN_INT *modified, FLOAT *threshold, FLOAT *ip, FLOAT *qtot)
{
  xc_gga_lb_set_params((xc_gga_type *)(*p), *modified, *threshold, *ip, *qtot);
}

void FC_FUNC_(xc_f90_gga_lb_modified, XC_F90_GGA_LB_MODIFIED)
     (void **p, FLOAT *rho, FLOAT *grho, FLOAT *r, FLOAT *dedd)
{
  xc_gga_lb_modified((xc_gga_type *)(*p), rho, grho, *r, dedd);
}


/* Hybrid GGAs */

void FC_FUNC_(xc_f90_hyb_gga_init, XC_F90_HYB_GGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_hyb_gga_type *hyb_gga_p;

  *p = malloc(sizeof(xc_hyb_gga_type));
  hyb_gga_p = (xc_hyb_gga_type *)(*p);
  xc_hyb_gga_init(hyb_gga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(hyb_gga_p->info);
}

void FC_FUNC_(xc_f90_hyb_gga_end, XC_F90_HYB_GGA_END)
     (void **p)
{
  xc_hyb_gga_end((xc_hyb_gga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_f90_hyb_gga, XC_F90_HYB_GGA)
     (void **p, FLOAT *rho, FLOAT *grho, 
      FLOAT *e, FLOAT *dedd, FLOAT *dedgd)
{
  xc_hyb_gga((xc_hyb_gga_type *)(*p), rho, grho, e, dedd, dedgd);
}

void FC_FUNC_(xc_f90_hyb_gga_exx_coef, XC_F90_HYB_GGA_EXX_COEF)
  (void **p, FLOAT *coef)
{
  *coef = xc_hyb_gga_exx_coef((xc_hyb_gga_type *)(*p));
}


/* meta-GGAs */

void FC_FUNC_(xc_f90_mgga_init, XC_F90_MGGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_mgga_type *mgga_p;

  *p = malloc(sizeof(xc_mgga_type));
  mgga_p = (xc_mgga_type *)(*p);
  xc_mgga_init(mgga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(mgga_p->info);
}

void FC_FUNC_(xc_f90_mgga_end, XC_F90_MGGA_END)
     (void **p)
{
  xc_mgga_end((xc_mgga_type *)(*p));
  free(*p);
}

void FC_FUNC_(xc_f90_mgga, XC_F90_MGGA)
  (void **p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
   FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau)
{
  xc_mgga((xc_mgga_type *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
}


/* LCAs */

void FC_FUNC_(xc_f90_lca_init, XC_F90_LCA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  xc_lca_type *lca_p;

  *p = malloc(sizeof(xc_lca_type));
  lca_p = (xc_lca_type *)(*p);
  xc_lca_init(lca_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(lca_p->info);
}

void FC_FUNC_(xc_f90_lca_end, XC_F90_LCA_END)
     (void **p)
{
  free(*p);
}

void FC_FUNC_(xc_f90_lca, XC_F90_LCA)
     (void **p, FLOAT *rho, FLOAT *v, 
      FLOAT *e, FLOAT *dedd, FLOAT *dedv)
{
  xc_lca((xc_lca_type *)(*p), rho, v, e, dedd, dedv);
}

#endif
