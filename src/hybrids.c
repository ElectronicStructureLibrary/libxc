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
xc_hyb_init_hybrid(xc_func_type *p, double alpha)
{
  p->hyb_number_terms = 1;
  
  p->hyb_type[0] = XC_HYB_FOCK;
  p->hyb_params[0][0] = alpha;
}

/* short-range hybrid. There are two parameters, the coefficient in
   from of the short-range Foch term, and the screening parameter
   omega */
void
xc_hyb_init_sr(xc_func_type *p, double beta, double omega)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_ERF_SR;
  p->hyb_params[0][0] = omega; /* for convenience omega is always the first */
  p->hyb_params[0][1] = beta;
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
  
  p->hyb_type[0] = XC_HYB_ERF_SR;
  p->hyb_params[0][0] = omega; /* for convenience omega is always the first */
  p->hyb_params[0][1] = beta;

  p->hyb_type[1] = XC_HYB_FOCK;
  p->hyb_params[1][0] = alpha;
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

/* van der Waals correction according to Dion2004_246401. One
   parameter only: Zab */
void
xc_hyb_init_vdw_df(xc_func_type *p, double Zab)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_VDW_DF;
  p->hyb_params[0][0] = Zab;
}

/* van der Waals correction according to Vydrov2010_244103. Two parameters, b and C */
void
xc_hyb_init_vdw_vv10(xc_func_type *p, double b, double C)
{
  p->hyb_number_terms = 1;

  p->hyb_type[0] = XC_HYB_VDW_VV10;
  p->hyb_params[0][0] = b;
  p->hyb_params[0][1] = C;
}

void
xc_hyb_init_vdw_df(xc_func_type *p, double Zab)
{
  int    hyb_type[1]  = {XC_HYB_VDW_DF};
  double hyb_omega[1] = {0.0};
  double hyb_coeff[1] = {Zab};

  xc_hyb_init(p, 1, hyb_type, hyb_coeff, hyb_omega);
}

void
xc_hyb_init_vdw_df(xc_func_type *p, double Zab)
{
  int    hyb_type[1]  = {XC_HYB_VDW_DF};
  double hyb_omega[1] = {0.0};
  double hyb_coeff[1] = {Zab};

  xc_hyb_init(p, 1, hyb_type, hyb_coeff, hyb_omega);
}

/* checks and returns the type of hybrid function */
int
xc_hyb_type(const xc_func_type *p)
{
  /* we check several specific types */
  if(p->hyb_number_terms == 0)
    return XC_HYB_NONE;

  if(p->hyb_number_terms == 1){
    /* normal hybrid */
    if(p->hyb_type[0] == XC_HYB_FOCK)
      return XC_HYB_HYBRID;

    /* range-separated hybrids */
    if(p->hyb_type[0] == XC_HYB_ERF_SR)
      return XC_HYB_CAM;
    if(p->hyb_type[0] == XC_HYB_YUKAWA_SR)
      return XC_HYB_CAMY;
    if(p->hyb_type[0] == XC_HYB_GAUSSIAN_SR)
      return XC_HYB_CAMG;

    /* van der Waals functional */
    if(p->hyb_type[0] == XC_HYB_VDW_DF)
      return XC_HYB_VDW;
  }

  if(p->hyb_number_terms == 2) {
    if(p->hyb_type[0] == XC_HYB_ERF_SR      && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAM;
    if(p->hyb_type[0] == XC_HYB_YUKAWA_SR   && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAMY;
    if(p->hyb_type[0] == XC_HYB_GAUSSIAN_SR && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAMG;
    
    if(p->hyb_type[0] == XC_HYB_PT2         && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_DOUBLE_HYBRID;
  }

  return XC_HYB_MIXTURE;
}


/*------------------------------------------------------*/
/* returns the mixing coefficient for the hybrid functions */
double
xc_hyb_exx_coef(const xc_func_type *p)
{
  assert(p!=NULL);
  assert(xc_hyb_type(p) == XC_HYB_HYBRID);

  return p->hyb_params[0][0];
}

/* returns the CAM parameters for (screened) hybrids */
void
xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta)
{
  assert(p != NULL);

  switch(xc_hyb_type(p)){
  case XC_HYB_HYBRID:
    *omega = 0.0;
    *beta  = 0.0;
    *alpha = p->hyb_params[0][0];
    break;
  case XC_HYB_CAM:
  case XC_HYB_CAMY:
  case XC_HYB_CAMG:
    if(p->hyb_number_terms == 1){ /* Short-range only hybrid */
      *omega = p->hyb_params[0][0];
      *beta  = p->hyb_params[0][1];
      *alpha = 0.0;
    }else{                        /* includes both short-range and normal hybrid */
      *omega = p->hyb_params[0][0];
      *beta  = p->hyb_params[0][1];
      *alpha = p->hyb_params[1][0];
    }
    break;
  default:
    fprintf(stderr, "Error in xc_hyb_cam_coef: unknown hybrid type");
    exit(1);
  }
}

/* returns the NLC parameters */
void
xc_hyb_vdw_vv10_coef(const xc_func_type *p, double *b, double *C)
{
  assert(p!=NULL);
  assert(p->hyb_type[0] == XC_HYB_VDW_VV10);
  
  *b = p->hyb_params[0][0];
  *C = p->hyb_params[0][1];;
}
