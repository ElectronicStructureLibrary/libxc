/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "config.h"

#ifdef HAVE_FORTRAN

#include "xc.h"
#include "string_f.h"

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

CC_FORTRAN_INT  XC_FC_FUNC(f90_info_provides, F90_INFO_PROVIDES)
     (void **info)
{
  return (CC_FORTRAN_INT) ((XC(func_info_type) *)(*info))->provides;
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
void XC_FC_FUNC(f90_lda_init, F90_LDA_INIT)
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
     (void **p, int *np, FLOAT *rho, 
      FLOAT *zk, FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, zk, vrho, v2rho2, v3rho3);
}

void XC_FC_FUNC(f90_lda_exc, F90_LDA_EXC)
     (void **p, int *np, FLOAT *rho,
      FLOAT *zk)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, zk, NULL, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_exc_vxc, F90_LDA_EXC_VXC)
     (void **p, int *np, FLOAT *rho, 
      FLOAT *zk, FLOAT *vrho)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, zk, vrho, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_vxc, F90_LDA_VXC)
     (void **p, int *np, FLOAT *rho, 
      FLOAT *vrho)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, NULL, vrho, NULL, NULL);
}

void XC_FC_FUNC(f90_lda_fxc, F90_LDA_FXC)
     (void **p, int *np, FLOAT *rho,
      FLOAT *v2rho2)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, NULL, NULL, v2rho2, NULL);
}

void XC_FC_FUNC(f90_lda_kxc, F90_LDA_KXC)
     (void **p, int *np, FLOAT *rho,
      FLOAT *v3rho3)
{
  XC(lda)((XC(lda_type) *)(*p), *np, rho, NULL, NULL, NULL, v3rho3);
}


/* Now come some special initializations */

/* parameter of Xalpha */
void XC_FC_FUNC(f90_lda_c_xalpha_set_par, F90_LDA_C_1D_CSC_SET_PAR)
  (void **p, FLOAT *alpha)
{
  XC(lda_c_xalpha_set_params)((XC(lda_type) *)(*p), *alpha);
}

/* parameter of CSC */
void XC_FC_FUNC(f90_lda_c_1d_csc_set_par, F90_LDA_C_1D_CSC_SET_PAR)
  (void **p, FLOAT *bb)
{
  XC(lda_c_1d_csc_set_params)((XC(lda_type) *)(*p), *bb);
}

/* parameter of PRM */
void XC_FC_FUNC(f90_lda_c_2d_prm_set_par, F90_LDA_C_2D_PRM_SET_PAR)
  (void **p, FLOAT *N)
{
  XC(lda_c_2d_prm_set_params)((XC(lda_type) *)(*p), *N);
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

void XC_FC_FUNC(f90_gga_exc_vxc, F90_GGA_EXC_VXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_exc_vxc)((XC(gga_type) *)(*p), rho, sigma, zk, vrho, vsigma);
}

void XC_FC_FUNC(f90_gga_vxc, F90_GGA_VXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(gga_vxc)((XC(gga_type) *)(*p), rho, sigma, vrho, vsigma);
}

void XC_FC_FUNC(f90_gga_fxc, F90_GGA_FXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(gga_fxc)((XC(gga_type) *)(*p), rho, sigma, v2rho2, v2rhosigma, v2sigma2);
}

/* the van Leeuwen & Baerends functional is special */
void XC_FC_FUNC(f90_gga_lb_set_par, F90_GGA_LB_SET_PAR)
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
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(hyb_gga)((XC(hyb_gga_type) *)(*p), rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

void XC_FC_FUNC(f90_hyb_gga_exc, F90_HYB_GGA_EXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk)
{
  XC(hyb_gga_exc)((XC(hyb_gga_type) *)(*p), rho, sigma, zk);
}

void XC_FC_FUNC(f90_hyb_gga_exc_vxc, F90_HYB_GGA_EXC_VXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma)
{
  XC(hyb_gga_exc_vxc)((XC(hyb_gga_type) *)(*p), rho, sigma, zk, vrho, vsigma);
}

void XC_FC_FUNC(f90_hyb_gga_vxc, F90_HYB_GGA_VXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *vrho, FLOAT *vsigma)
{
  XC(hyb_gga_vxc)((XC(hyb_gga_type) *)(*p), rho, sigma, vrho, vsigma);
}

void XC_FC_FUNC(f90_hyb_gga_fxc, F90_HYB_GGA_FXC)
     (void **p, FLOAT *rho, FLOAT *sigma, 
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  XC(hyb_gga_fxc)((XC(hyb_gga_type) *)(*p), rho, sigma, v2rho2, v2rhosigma, v2sigma2);
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
     (void **p, FLOAT *rho, FLOAT *sigma, FLOAT *lapl_rho, FLOAT *tau,
      FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau,
      FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga)((XC(mgga_type) *)(*p), rho, sigma, lapl_rho, tau, 
	   zk, vrho, vsigma, vlapl_rho, vtau,
	   v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}

void XC_FC_FUNC(f90_mgga_exc, F90_MGGA_EXC)
     (void **p, FLOAT *rho, FLOAT *sigma, FLOAT *lapl_rho, FLOAT *tau, 
      FLOAT *zk)
{
  XC(mgga_exc)((XC(mgga_type) *)(*p), rho, sigma, lapl_rho, tau, zk);
}

void XC_FC_FUNC(f90_mgga_exc_vxc, F90_MGGA_EXC_VXC)
  (void **p, FLOAT *rho, FLOAT *sigma, FLOAT *lapl_rho, FLOAT *tau,
   FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga_exc_vxc)((XC(mgga_type) *)(*p), rho, sigma, lapl_rho, tau, zk, vrho, vsigma, vlapl_rho, vtau);
}

void XC_FC_FUNC(f90_mgga_vxc, F90_MGGA_VXC)
  (void **p, FLOAT *rho, FLOAT *sigma, FLOAT *lapl_rho, FLOAT *tau,
   FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho, FLOAT *vtau)
{
  XC(mgga_vxc)((XC(mgga_type) *)(*p), rho, sigma, lapl_rho, tau, vrho, vsigma, vlapl_rho, vtau);
}

void XC_FC_FUNC(f90_mgga_fxc, F90_MGGA_FXC)
  (void **p, FLOAT *rho, FLOAT *sigma, FLOAT *lapl_rho, FLOAT *tau,
   FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  XC(mgga_fxc)((XC(mgga_type) *)(*p), rho, sigma, lapl_rho, tau, v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
}

/* parameter of TP09 */
void XC_FC_FUNC(f90_mgga_x_tb09_set_par, F90_MGGA_X_TB09_SET_PAR)
  (void **p, FLOAT *cc)
{
  XC(mgga_x_tb09_set_params)((XC(mgga_type) *)(*p), *cc);
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
