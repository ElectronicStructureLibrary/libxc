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
}

/* Double precision interfaces */

void FC_FUNC_(xc_f90_lda_dp, XC_F90_LDA_DP)
  (void **p, double *rho, double *exc, double *vxc, double *fxc, double *kxc)
{
  xc_lda((xc_lda_type *)(*p), rho, exc, vxc, fxc, kxc);
}

void FC_FUNC_(xc_f90_lda_vxc_dp, XC_F90_LDA_VXC_DP)
     (void **p, double *rho, double *e, double *v)
{
  xc_lda_vxc((xc_lda_type *)(*p), rho, e, v);
}

void FC_FUNC_(xc_f90_lda_fxc_dp, XC_F90_LDA_FXC_DP)
     (void **p, double *rho, double *fxc)
{
  xc_lda_fxc((xc_lda_type *)(*p), rho, fxc);
}

void FC_FUNC_(xc_f90_lda_kxc_dp, XC_F90_LDA_KXC_DP)
     (void **p, double *rho, double *kxc)
{
  xc_lda_kxc((xc_lda_type *)(*p), rho, kxc);
}

/* Single precision interfaces */

void FC_FUNC_(xc_f90_lda_sp, XC_F90_LDA_SP)
  (void **p, float *rho, float *exc, float *vxc, float *fxc, float *kxc)
{
  xc_lda_sp((xc_lda_type *)(*p), rho, exc, vxc, fxc, kxc);
}

void FC_FUNC_(xc_f90_lda_vxc_sp, XC_F90_LDA_VXC_SP)
     (void **p, float *rho, float *e, float *v)
{
  xc_lda_vxc_sp((xc_lda_type *)(*p), rho, e, v);
}

void FC_FUNC_(xc_f90_lda_fxc_sp, XC_F90_LDA_FXC_SP)
     (void **p, float *rho, float *fxc)
{
  xc_lda_fxc_sp((xc_lda_type *)(*p), rho, fxc);
}

void FC_FUNC_(xc_f90_lda_kxc_sp, XC_F90_LDA_KXC_SP)
     (void **p, float *rho, float *kxc)
{
  xc_lda_kxc_sp((xc_lda_type *)(*p), rho, kxc);
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
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, double *alpha)
{
  xc_lda_type *lda_p;

  assert((int) (*functional) == XC_LDA_C_XALPHA);

  *p = malloc(sizeof(xc_lda_type));
  lda_p = (xc_lda_type *)(*p);
  xc_lda_c_xalpha_init(lda_p, (int) (*nspin), (int) (*dim), *alpha);
  *info = (void *)(lda_p->info);
}

/* single precision version of the previous */
void FC_FUNC_(xc_f90_lda_c_xalpha_init_sp, XC_F90_LDA_C_XALPHA_INIT_SP)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, float *alpha)
{
  double dalpha;
  dalpha = alpha[0];
  FC_FUNC_(xc_f90_lda_c_xalpha_init, XC_F90_LDA_C_XALPHA_INIT)
    (p, info, functional, nspin, dim, &dalpha);
}


/* GGAs */

void FC_FUNC_(xc_f90_gga_init_, XC_F90_GGA_INIT_)
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
}

void FC_FUNC_(xc_f90_gga_dp, XC_F90_GGA_DP)
     (void **p, double *rho, double *grho, 
      double *e, double *dedd, double *dedgd)
{
  xc_gga((xc_gga_type *)(*p), rho, grho, e, dedd, dedgd);
}

void FC_FUNC_(xc_f90_gga_sp, XC_F90_GGA_SP)
     (void **p, float *rho, float *grho, 
      float *e, float *dedd, float *dedgd)
{
  xc_gga_sp((xc_gga_type *)(*p), rho, grho, e, dedd, dedgd);
}


/* the van Leeuwen & Baerends functional is special */
void FC_FUNC_(xc_f90_gga_lb_init, XC_F90_GGA_LB_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *modified, double *threshold)
{
  xc_gga_type *gga_p;

  assert(*functional == XC_GGA_XC_LB);

  *p = malloc(sizeof(xc_gga_type));
  gga_p = (xc_gga_type *)(*p);
  xc_gga_lb_init(gga_p, (int) (*nspin), (int) (*modified), (int) (*threshold));
  *info = (void *)(gga_p->info);
}

void FC_FUNC_(xc_f90_gga_lb_init_sp, XC_F90_GGA_LB_INIT_SP)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *modified, float *threshold)
{
  double dthreshold = threshold[0];
  FC_FUNC_(xc_f90_gga_lb_init, XC_F90_GGA_LB_INIT)
    (p, info, functional, nspin, modified, &dthreshold);
}

void FC_FUNC_(xc_f90_gga_lb_dp, XC_F90_GGA_LB_DP)
     (void **p, double *rho, double *grho, double *r, double *ip, double *qtot,
      double *dedd)
{
  xc_gga_lb((xc_gga_type *)(*p), rho, grho, *r, *ip, *qtot, dedd);
}

void FC_FUNC_(xc_f90_gga_lb_sp, XC_F90_GGA_LB_SP)
     (void **p, float *rho, float *grho, float *r, float *ip, float *qtot,
      float *dedd)
{
  xc_gga_lb_sp((xc_gga_type *)(*p), rho, grho, *r, *ip, *qtot, dedd);
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

void FC_FUNC_(xc_f90_mgga_dp, XC_F90_MGGA_DP)
  (void **p, double *rho, double *grho, double *tau,
   double *e, double *dedd, double *dedgd, double *dedtau)
{
  xc_mgga((xc_mgga_type *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
}

void FC_FUNC_(xc_f90_mgga_sp, XC_F90_MGGA_SP)
  (void **p, float *rho, float *grho, float *tau,
   float *e, float *dedd, float *dedgd, float *dedtau)
{
  xc_mgga_sp((xc_mgga_type *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
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

void FC_FUNC_(xc_f90_lca_dp, XC_F90_LCA_DP)
     (void **p, double *rho, double *v, 
      double *e, double *dedd, double *dedv)
{
  xc_lca((xc_lca_type *)(*p), rho, v, e, dedd, dedv);
}

void FC_FUNC_(xc_f90_lca_sp, XC_F90_LCA_SP)
  (void **p, float *rho, float *v, 
   float *e, float *dedd, float *dedv)
{
  xc_lca_sp((xc_lca_type *)(*p), rho, v, e, dedd, dedv);
}

#endif
