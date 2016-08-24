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
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "xc.h"
#include "funcs_key.c"

extern XC(func_info_type) 
  *XC(lda_known_funct)[], 
  *XC(gga_known_funct)[],
  *XC(hyb_gga_known_funct)[],
  *XC(mgga_known_funct)[],
  *XC(hyb_mgga_known_funct)[];


/*------------------------------------------------------*/
int XC(functional_get_number)(const char *name)
{
  int ii;
  int key=-1;
  const char *p;

  /* Does name begin with xc_? */
  if(strncasecmp(name,"XC_",3) == 0) {
    p=name+3;
  } else {
    p=name;
  }

  for(ii=0;;ii++){
    if(XC(functional_keys)[ii].number == -1)
      break;
    if(strcasecmp(XC(functional_keys)[ii].name, p) == 0){
      key = XC(functional_keys)[ii].number;
      break;
    }
  }
  
  return key;
}


/*------------------------------------------------------*/
char *XC(functional_get_name)(const int number)
{
  int ii;

  for(ii=0;;ii++){
    if(XC(functional_keys)[ii].number == -1)
      return NULL;
    if(XC(functional_keys)[ii].number == number)
      /* return duplicated: caller has the responsibility to dealloc string */
      return strdup(XC(functional_keys)[ii].name);
  }
}


/*------------------------------------------------------*/
int XC(family_from_id)(int id, int *family, int *number)
{
  int ii;

  /* first let us check if it is an LDA */
  for(ii=0; XC(lda_known_funct)[ii]!=NULL; ii++){
    if(XC(lda_known_funct)[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_LDA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_LDA;
    }
  }

  /* or is it a GGA? */
  for(ii=0; XC(gga_known_funct)[ii]!=NULL; ii++){
    if(XC(gga_known_funct)[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_GGA;
    }
  }

  /* or is it a hybrid GGA? */
  for(ii=0; XC(hyb_gga_known_funct)[ii]!=NULL; ii++){
    if(XC(hyb_gga_known_funct)[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_GGA;
    }
  }

  /* or is it a meta GGA? */
  for(ii=0; XC(mgga_known_funct)[ii]!=NULL; ii++){
    if(XC(mgga_known_funct)[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_MGGA;
    }
  }

  /* or is it a hybrid meta GGA? */
  for(ii=0; XC(hyb_mgga_known_funct)[ii]!=NULL; ii++){
    if(XC(hyb_mgga_known_funct)[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_MGGA;
    }
  }

  return XC_FAMILY_UNKNOWN;
}

/*------------------------------------------------------*/
XC(func_type) *XC(func_alloc)()
{
  XC(func_type) *func;

  func = (XC(func_type) *) malloc (sizeof (XC(func_type)));
  return func;
}

/*------------------------------------------------------*/
int XC(func_init)(XC(func_type) *func, int functional, int nspin)
{
  int number;

  assert(func != NULL);
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);

  /* initialize structure */
  func->nspin      = nspin;
  func->params     = NULL;
  func->func       = 0;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix_coef   = NULL;
  func->cam_omega = func->cam_alpha = func->cam_beta = 0.0;
  func->nlc_b = func->nlc_C = 0.0;

  switch(XC(family_from_id)(functional, NULL, &number)){
  case(XC_FAMILY_LDA):
    func->info = XC(lda_known_funct)[number];
    break;

  case(XC_FAMILY_GGA):
    func->info = XC(gga_known_funct)[number];
    break;

  case(XC_FAMILY_HYB_GGA):
    func->info = XC(hyb_gga_known_funct)[number];
    break;

  case(XC_FAMILY_MGGA):
    func->info = XC(mgga_known_funct)[number];
    break;

  case(XC_FAMILY_HYB_MGGA):
    func->info = XC(hyb_mgga_known_funct)[number];
    break;

  default:
    return -2; /* family not found */
  }

  /* setup internal counters */
  switch(XC(family_from_id)(functional, NULL, &number)){
  case(XC_FAMILY_MGGA):
  case(XC_FAMILY_HYB_MGGA):
    func->n_tau  = func->n_vtau = func->nspin;
    func->n_lapl = func->n_vlapl = func->nspin;
    if(func->nspin == XC_UNPOLARIZED){
      func->n_v2tau2 = func->n_v2lapl2 = 1;
      func->n_v2rhotau = func->n_v2rholapl = func->n_v2lapltau = 1;
      func->n_v2sigmatau = func->n_v2sigmalapl = 1;
    }else{
      func->n_v2tau2 = func->n_v2lapl2 = 3;
      func->n_v2rhotau = func->n_v2rholapl = func->n_v2lapltau = 4;
      func->n_v2sigmatau = func->n_v2sigmalapl = 6;
    }

  case(XC_FAMILY_GGA):
  case(XC_FAMILY_HYB_GGA):
    if(func->nspin == XC_UNPOLARIZED){
      func->n_sigma  = func->n_vsigma = 1;
      func->n_v2rhosigma  = func->n_v2sigma2 = 1;
      func->n_v3rho2sigma = func->n_v3rhosigma2 = func->n_v3sigma3 = 1;
    }else{
      func->n_sigma      = func->n_vsigma = 3;
      func->n_v2rhosigma = func->n_v2sigma2 = 6;

      func->n_v3rho2sigma = 9;
      func->n_v3rhosigma2 = 12;
      func->n_v3sigma3    = 10;
    }

  case(XC_FAMILY_LDA):
    func->n_rho = func->n_vrho = func->nspin;
    func->n_zk  = 1;
    if(func->nspin == XC_UNPOLARIZED){
      func->n_v2rho2 = func->n_v3rho3 = 1;
    }else{
      func->n_v2rho2 = 3;
      func->n_v3rho3 = 4;
    }
  }

  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);

  /* see if we need to initialize the external parameters */
  if(func->info->n_ext_params > 0)
    func->info->set_ext_params(func, NULL);

  return 0;
}


/*------------------------------------------------------*/
void XC(func_end)(XC(func_type) *func)
{
  assert(func != NULL && func->info != NULL);

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(func);

  /* terminate any auxiliary functional */
  if(func->n_func_aux > 0){
    int ii;

    for(ii=0; ii<func->n_func_aux; ii++){
      XC(func_end)(func->func_aux[ii]);
      free(func->func_aux[ii]);
    }
    free(func->func_aux);
    func->n_func_aux = 0;
  }

  if(func->mix_coef != NULL){
    free(func->mix_coef);
    func->mix_coef = NULL;
  }

  /* deallocate any used parameter */
  if(func->params != NULL){
    free(func->params);
    func->params = NULL;
  }

  func->info = NULL;  
}

/*------------------------------------------------------*/
void  XC(func_free)(XC(func_type) *p)
{
  free(p);
}

/*------------------------------------------------------*/
const XC(func_info_type) *XC(func_get_info)(const XC(func_type) *p)
{
  return p->info;
}

/*------------------------------------------------------*/
void XC(func_set_ext_params)(XC(func_type) *p, double *ext_params)
{
  assert(p->info->n_ext_params > 0);
  p->info->set_ext_params(p, ext_params);
}

/* returns the mixing coefficient for the hybrid GGAs */
FLOAT XC(hyb_exx_coef)(const XC(func_type) *p)
{
  assert(p!=NULL);
 
  return p->cam_alpha;
}

/* returns the CAM parameters for screened hybrids */
void XC(hyb_cam_coef)(const XC(func_type) *p, FLOAT *omega, FLOAT *alpha, FLOAT *beta)
{
  assert(p!=NULL);

  *omega = p->cam_omega;
  *alpha = p->cam_alpha;
  *beta  = p->cam_beta;
}

/* returns the NLC parameters */
void XC(nlc_coef)(const XC(func_type) *p, FLOAT *nlc_b, FLOAT *nlc_C)
{
  assert(p!=NULL);

  *nlc_b = p->nlc_b;
  *nlc_C = p->nlc_C;
}
