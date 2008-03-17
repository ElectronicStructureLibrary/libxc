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

CC_FORTRAN_INT XC_FC_FUNC(f90_info_number, F90_INFO_NUMBER)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->number;
}


CC_FORTRAN_INT XC_FC_FUNC(f90_info_kind, F90_INFO_KIND)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->kind;
}


void XC_FC_FUNC(f90_info_name, F90_INFO_NAME)
     (void **info, STR_F_TYPE s STR_ARG1)
{
  TO_F_STR1(((XC(func_info_type) *)(*info))->name, s);
}


CC_FORTRAN_INT  XC_FC_FUNC(f90_info_family, F90_INFO_FAMILY)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->family;
}

void XC_FC_FUNC(f90_info_ref, F90_INFO_REF)
     (void **info, char **s, STR_F_TYPE ref_f STR_ARG1)
{
  char *c, ref[256]; /* hopefully no ref is longer than 256 characters ;) */
  XC(func_info_type) *func_p = (XC(func_info_type) *)(*info);

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
CC_FORTRAN_INT  XC_FC_FUNC(f90_family_from_id, F90_FAMILY_FROM_ID)
  (CC_FORTRAN_INT  *functional)
{
  return (CC_FORTRAN_INT) XC(family_from_id)((int) (*functional));
}


/* LDAs */

/* Standard initialization */
void XC_FC_FUNC(f90_lda_init_, F90_LDA_INIT_)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(lda_type) *lda_p;
  
  *p = malloc(sizeof(XC(lda_type)));
  lda_p = (XC(lda_type) *)(*p);
  XC(lda_init)(lda_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(lda_p->info);
}

void XC_FC_FUNC(f90_lda_end, F90_LDA_END)
     (void **p)
{
  /* xc_lda_end does not exist */
  free(*p);
  *p = NULL;
}

void XC_FC_FUNC(f90_lda, F90_LDA)
  (void **p, FLOAT *rho, 
   FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
{
  XC(lda)((XC(lda_type) *)(*p), rho, zk, vrho, v2rho2, v3rho3);
}

void XC_FC_FUNC(f90_lda_exc, F90_LDA_EXC)
     (void **p, FLOAT *rho, FLOAT *zk)
{
  XC(lda_exc)((XC(lda_type) *)(*p), rho, zk);
}

void XC_FC_FUNC(f90_lda_vxc, F90_LDA_VXC)
     (void **p, FLOAT *rho, FLOAT *zk, FLOAT *vrho)
{
  XC(lda_vxc)((XC(lda_type) *)(*p), rho, zk, vrho);
}

void XC_FC_FUNC(f90_lda_fxc, F90_LDA_FXC)
     (void **p, FLOAT *rho, FLOAT *v2rho2)
{
  XC(lda_fxc)((XC(lda_type) *)(*p), rho, v2rho2);
}

void XC_FC_FUNC(f90_lda_kxc, F90_LDA_KXC)
     (void **p, FLOAT *rho, FLOAT *v3rho3)
{
  XC(lda_kxc)((XC(lda_type) *)(*p), rho, v3rho3);
}


/* Now come some special initializations */

/* exchange in the LDA */
void XC_FC_FUNC(f90_lda_x_init, F90_LDA_X_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, CC_FORTRAN_INT *irel)
{
  XC(lda_type) *lda_p;

  assert(*functional == XC_LDA_X);

  *p = malloc(sizeof(XC(lda_type)));
  lda_p = (XC(lda_type) *)(*p);
  XC(lda_x_init)(lda_p, (int) (*nspin), (int) (*dim), (int) (*irel));
  *info = (void *)(lda_p->info);
}

/* Slater's Xalpha */
void XC_FC_FUNC(f90_lda_c_xalpha_init, F90_LDA_C_XALPHA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, 
      CC_FORTRAN_INT *nspin, CC_FORTRAN_INT *dim, FLOAT *alpha)
{
  XC(lda_type) *lda_p;

  assert((int) (*functional) == XC_LDA_C_XALPHA);

  *p = malloc(sizeof(XC(lda_type)));
  lda_p = (XC(lda_type) *)(*p);
  XC(lda_c_xalpha_init)(lda_p, (int) (*nspin), (int) (*dim), *alpha);
  *info = (void *)(lda_p->info);
}


/* GGAs */

void XC_FC_FUNC(f90_gga_init, F90_GGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(gga_type) *gga_p;

  *p = malloc(sizeof(XC(gga_type)));
  gga_p = (XC(gga_type) *)(*p);
  XC(gga_init)(gga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(gga_p->info);
}

void XC_FC_FUNC(f90_gga_end, F90_GGA_END)
     (void **p)
{
  XC(gga_end)((XC(gga_type) *)(*p));
  free(*p);
  *p = NULL;
}

void XC_FC_FUNC(f90_gga, F90_GGA)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga)((XC(gga_type) *)(*p), rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

void XC_FC_FUNC(f90_gga_exc, F90_GGA_EXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk)
{
  XC(gga_exc)((XC(gga_type) *)(*p), rho, sigma, zk);
}

void XC_FC_FUNC(f90_gga_vxc, F90_GGA_VXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_vxc)((XC(gga_type) *)(*p), rho, sigma, zk, vrho, vsigma);
}

void XC_FC_FUNC(f90_gga_fxc, F90_GGA_FXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga_fxc)((XC(gga_type) *)(*p), rho, sigma, v2rho2, v2rhosigma, v2sigma2);
}

/* the van Leeuwen & Baerends functional is special */
void XC_FC_FUNC(f90_gga_lb_set_params, F90_GGA_LB_SET_PARAMS)
  (void **p, CC_FORTRAN_INT *modified, FLOAT *threshold, FLOAT *ip, FLOAT *qtot)
{
  XC(gga_lb_set_params)((XC(gga_type) *)(*p), *modified, *threshold, *ip, *qtot);
}

void XC_FC_FUNC(f90_gga_lb_modified, F90_GGA_LB_MODIFIED)
     (void **p, FLOAT *rho, FLOAT *grho, FLOAT *r, FLOAT *dedd)
{
  XC(gga_lb_modified)((XC(gga_type) *)(*p), rho, grho, *r, dedd);
}


/* Hybrid GGAs */

void XC_FC_FUNC(f90_hyb_gga_init, F90_HYB_GGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(hyb_gga_type) *hyb_gga_p;

  *p = malloc(sizeof(XC(hyb_gga_type)));
  hyb_gga_p = (XC(hyb_gga_type) *)(*p);
  XC(hyb_gga_init)(hyb_gga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(hyb_gga_p->info);
}

void XC_FC_FUNC(f90_hyb_gga_end, F90_HYB_GGA_END)
     (void **p)
{
  XC(hyb_gga_end)((XC(hyb_gga_type) *)(*p));
  free(*p);
}

void XC_FC_FUNC(f90_hyb_gga, F90_HYB_GGA)
     (void **p, FLOAT *rho, FLOAT *grho, 
      FLOAT *e, FLOAT *dedd, FLOAT *dedgd)
{
  XC(hyb_gga)((XC(hyb_gga_type) *)(*p), rho, grho, e, dedd, dedgd);
}

void XC_FC_FUNC(f90_hyb_gga_exx_coef, F90_HYB_GGA_EXX_COEF)
  (void **p, FLOAT *coef)
{
  *coef = XC(hyb_gga_exx_coef)((XC(hyb_gga_type) *)(*p));
}


/* meta-GGAs */

void XC_FC_FUNC(f90_mgga_init, F90_MGGA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(mgga_type) *mgga_p;

  *p = malloc(sizeof(XC(mgga_type)));
  mgga_p = (XC(mgga_type) *)(*p);
  XC(mgga_init)(mgga_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(mgga_p->info);
}

void XC_FC_FUNC(f90_mgga_end, F90_MGGA_END)
     (void **p)
{
  XC(mgga_end)((XC(mgga_type) *)(*p));
  free(*p);
}

void XC_FC_FUNC(f90_mgga, F90_MGGA)
  (void **p, FLOAT *rho, FLOAT *grho, FLOAT *tau,
   FLOAT *e, FLOAT *dedd, FLOAT *dedgd, FLOAT *dedtau)
{
  XC(mgga)((XC(mgga_type) *)(*p), rho, grho, tau, e, dedd, dedgd, dedtau);
}


/* LCAs */

void XC_FC_FUNC(f90_lca_init, F90_LCA_INIT)
     (void **p, void **info, CC_FORTRAN_INT *functional, CC_FORTRAN_INT *nspin)
{
  XC(lca_type) *lca_p;

  *p = malloc(sizeof(XC(lca_type)));
  lca_p = (XC(lca_type) *)(*p);
  XC(lca_init)(lca_p, (int) (*functional), (int) (*nspin));
  *info = (void *)(lca_p->info);
}

void XC_FC_FUNC(f90_lca_end, F90_LCA_END)
     (void **p)
{
  free(*p);
}

void XC_FC_FUNC(f90_lca, F90_LCA)
     (void **p, FLOAT *rho, FLOAT *v, 
      FLOAT *e, FLOAT *dedd, FLOAT *dedv)
{
  XC(lca)((XC(lca_type) *)(*p), rho, v, e, dedd, dedv);
}

#endif
