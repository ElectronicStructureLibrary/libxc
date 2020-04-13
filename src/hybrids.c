/*
 Copyright (C) 2020 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

void xc_hyb_init(xc_func_type *p, int n_terms, const int *type,
                 const double *coeff, const double *omega)
{
  int ii;
  
  p->hyb_number_terms = n_terms;

  p->hyb_type  = (int    *) libxc_malloc(n_terms*sizeof(int));
  p->hyb_coeff = (double *) libxc_malloc(n_terms*sizeof(double));
  p->hyb_omega = (double *) libxc_malloc(n_terms*sizeof(double));

  for(ii=0; ii<n_terms; ii++){
    p->hyb_type[ii]  = type[ii];
    p->hyb_coeff[ii] = coeff[ii];
    p->hyb_omega[ii] = omega[ii];
  }
}

/* some specializations */
void
xc_hyb_init_hybrid(xc_func_type *p, double alpha)
{
  int    hyb_type[1]  = {XC_HYB_FOCK};
  double hyb_omega[1] = {0.0};
  double hyb_coeff[1] = {alpha};

  xc_hyb_init(p, 1, hyb_type, hyb_coeff, hyb_omega);
}

void
xc_hyb_init_sr(xc_func_type *p, double omega, double beta)
{
  int    hyb_type[1]  = {XC_HYB_ERF_SR};
  double hyb_omega[1] = {omega};
  double hyb_coeff[1] = {beta};

  xc_hyb_init(p, 1, hyb_type, hyb_coeff, hyb_omega);
}

void
xc_hyb_init_cam(xc_func_type *p, double omega, double alpha, double beta)
{
  int hyb_type[2]     = {XC_HYB_ERF_SR, XC_HYB_FOCK};
  double hyb_omega[2] = {omega, 0.0};
  double hyb_coeff[2] = {beta, alpha};

  xc_hyb_init(p, 2, hyb_type, hyb_coeff, hyb_omega);
}

void
xc_hyb_init_camy(xc_func_type *p, double omega, double alpha, double beta)
{
  int hyb_type[2]     = {XC_HYB_YUKAWA_SR, XC_HYB_FOCK};
  double hyb_omega[2] = {omega, 0.0};
  double hyb_coeff[2] = {beta, alpha};

  xc_hyb_init(p, 2, hyb_type, hyb_coeff, hyb_omega);
}

void
xc_hyb_init_camg(xc_func_type *p, double omega, double alpha, double beta)
{
  int hyb_type[2]     = {XC_HYB_GAUSSIAN_SR, XC_HYB_FOCK};
  double hyb_omega[2] = {omega, 0.0};
  double hyb_coeff[2] = {beta, alpha};

  xc_hyb_init(p, 2, hyb_type, hyb_coeff, hyb_omega);
}

/* checks and returns the type of hybrid function */
int
xc_hyb_type(const xc_func_type *p)
{
  /* we check several specific types */
  if(p->hyb_number_terms == 0)
    return XC_HYB_NONE;

  if(p->hyb_number_terms == 1){
    /* some GGAs are screened and keep omega in this structure */
    if(p->hyb_type[0] == XC_HYB_NONE)
      return XC_HYB_SEMILOCAL;

    /* normal hybrid */
    if(p->hyb_type[0] == XC_HYB_FOCK)
      return XC_HYB_HYBRID;

    if(p->hyb_type[0] == XC_HYB_ERF_SR)
      return XC_HYB_SHORT_RANGE;
  }

  if(p->hyb_number_terms == 2){
    if(p->hyb_type[0] == XC_HYB_ERF_SR      && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAM;
    if(p->hyb_type[0] == XC_HYB_YUKAWA_SR   && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAMY;
    if(p->hyb_type[0] == XC_HYB_GAUSSIAN_SR && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_CAMG;
    if(p->hyb_type[0] == XC_HYB_PT2         && p->hyb_type[1] == XC_HYB_FOCK)
      return XC_HYB_DOUBLE;
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
  
  return p->hyb_coeff[0];
}

/* returns the CAM parameters for screened hybrids */
void
xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta)
{
  assert(p!=NULL);
  assert(xc_hyb_type(p) == XC_HYB_CAM || xc_hyb_type(p) == XC_HYB_CAMY || xc_hyb_type(p) == XC_HYB_CAMG);

  *omega = p->hyb_omega[0];
  *beta  = p->hyb_coeff[0];
  *alpha = p->hyb_coeff[1];
}
