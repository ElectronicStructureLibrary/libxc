/*
 Copyright (C) 2020 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* standard hybrid functional that is defined by a single parameter
   alpha, specifically the coefficient in front of the Fock term */
void
xc_hyb_init_fock(xc_func_type *p, double alpha)
{
  p->hyb_number_terms = 1;
  
  p->hyb_type[0] = XC_HYB_FOCK;
  p->hyb_params[0].fock.alpha = alpha;
}

/* short-range hybrid. There are two parameters, the coefficient in
   from of the short-range Foch term, and the screening parameter
   omega */
void
xc_hyb_init_sr(xc_func_type *p, double beta, double omega)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_ERF_SR;
  p->hyb_params[0].sr.beta  = beta;
  p->hyb_params[0].sr.omega = omega;
}

/* Coulomb attenuated hybrid, based on the erf attenuation function.
     N.B. Different conventions for alpha and beta can be found in
     literature. In the convention used in libxc, at short range the
     fraction of exact exchange is alpha + beta, while at long range
     it is alpha. */
void
xc_hyb_init_cam(xc_func_type *p, double alpha, double beta, double omega)
{
  p->hyb_number_terms = 2;
  
  p->hyb_type[0] = XC_HYB_FOCK;
  p->hyb_params[0].fock.alpha = alpha;

  p->hyb_type[1] = XC_HYB_ERF_SR;
  p->hyb_params[1].sr.beta  = beta;
  p->hyb_params[1].sr.omega = omega;
}

/* Coulomb attenuated hybrid, based on the yukawa attenuation
   function. */
void
xc_hyb_init_camy(xc_func_type *p, double alpha, double beta, double omega)
{
  xc_hyb_init_cam(p, alpha, beta, omega);
  
  p->hyb_type[0] = XC_HYB_YUKAWA_SR;
}

/* Coulomb attenuated hybrid, based on the gaussian attenuation
   function. */
void
xc_hyb_init_camg(xc_func_type *p, double alpha, double beta, double omega)
{
  xc_hyb_init_cam(p, alpha, beta, omega);
  
  p->hyb_type[0] = XC_HYB_GAUSSIAN_SR;
}

/* van der Waals correction according to Grimme */
void
xc_hyb_init_vdw_d(xc_func_type *p, int type, double s6, double alpha, double r0)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = type;
  p->hyb_params[0].d.s6    = s6;
  p->hyb_params[0].d.alpha = alpha;
  p->hyb_params[0].d.r0    = r0;
}

/* van der Waals correction according to Dion2004_246401. One
   parameter only: Zab */
void
xc_hyb_init_vdw_df(xc_func_type *p, double delta, double Zab)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_VDW_DF;
  p->hyb_params[0].df.delta = delta;
  p->hyb_params[0].df.Zab   = Zab;
}

/* van der Waals correction according to Vydrov2010_244103. Two
   parameters, b and C */
void
xc_hyb_init_vdw_vv10(xc_func_type *p, double delta, double b, double C)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_VDW_VV10;

  p->hyb_params[0].vv10.delta = delta;
  p->hyb_params[0].vv10.b     = b;
  p->hyb_params[0].vv10.C     = C;
}

/* checks and returns the type of hybrid function */
int
xc_hyb_type(const xc_func_type *p)
{
  int result = 0;
  for(int i=0; i < p->hyb_number_terms; i++) {
    result = result | p->hyb_type[i];
  }
  return result;
}


/*------------------------------------------------------*/
/* returns the mixing coefficient for the hybrid functions */
double
xc_hyb_exx_coef(const xc_func_type *p)
{
  double cx=0.0;
  assert(p!=NULL);

  return p->hyb_params[0].fock.alpha;
}

/* returns the CAM parameters for (screened) hybrids */
void
xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta)
{
  /* Flags to check that the functional doesn't have duplicate contributions */
  int nfock=0;
  int nrangesep=0;
  /* Check pointer */
  assert(p != NULL);
  /* Initialize */
  *alpha = 0.0;
  *beta = 0.0;
  *omega = 0.0;

  /* Collect the parameters */
  for(int i=0;i<p->hyb_number_terms;i++) {
    switch(p->hyb_type[i]){
    case XC_HYB_FOCK:
      *alpha = p->hyb_params[i].fock.alpha;
      nfock++;
      break;
    case XC_HYB_ERF_SR:
    case XC_HYB_YUKAWA_SR:
    case XC_HYB_GAUSSIAN_SR:
      *beta  = p->hyb_params[i].sr.beta;
      *omega = p->hyb_params[i].sr.omega;
      nrangesep++;
      break;
    }
  }

  /* Ensure the functional did not have multiple components of the
     same type (this case needs to be handled by specialized code in
     the calling program as xc_hyb_cam_coef is mainly for backwards
     compatibility) */
  if(nfock>1) {
    fprintf(stderr, "Error in xc_hyb_cam_coef: functional contains multiple Fock matrix contributions\n");
  }
  if(nrangesep>1) {
    fprintf(stderr, "Error in xc_hyb_cam_coef: functional contains multiple range-separated contributions\n");
  }
}

/* returns the vdw parameters */
void
xc_hyb_vdw_df_coef(const xc_func_type *p, double *delta, double *Zab)
{
  assert(p!=NULL);
  assert(p->hyb_type[0] == XC_HYB_VDW_DF);

  *delta = p->hyb_params[0].df.delta;
  *Zab   = p->hyb_params[0].df.Zab;
}

void
xc_hyb_vdw_vv10_coef(const xc_func_type *p, double *delta, double *b, double *C)
{
  assert(p!=NULL);
  assert(p->hyb_type[0] == XC_HYB_VDW_VV10);
  
  *delta = p->hyb_params[0].vv10.delta;
  *b     = p->hyb_params[0].vv10.b;
  *C     = p->hyb_params[0].vv10.C;
}
