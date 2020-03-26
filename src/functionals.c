/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "xc.h"
#include "funcs_key.c"
#include <string.h>
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

extern xc_func_info_type 
  *xc_lda_known_funct[],
  *xc_hyb_lda_known_funct[],
  *xc_gga_known_funct[],
  *xc_hyb_gga_known_funct[],
  *xc_mgga_known_funct[],
  *xc_hyb_mgga_known_funct[];


/*------------------------------------------------------*/
int xc_functional_get_number(const char *name)
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
    if(xc_functional_keys[ii].number == -1)
      break;
    if(strcasecmp(xc_functional_keys[ii].name, p) == 0){
      key = xc_functional_keys[ii].number;
      break;
    }
  }
  
  return key;
}


/*------------------------------------------------------*/
char *xc_functional_get_name(int number)
{
  int ii;
  char *p;

  for(ii=0;;ii++){
    if(xc_functional_keys[ii].number == -1)
      return NULL;
    if(xc_functional_keys[ii].number == number) {
      /* return duplicated: caller has the responsibility to dealloc string.
         Do this the old way since strdup and strndup aren't C standard. */
      p = (char *) libxc_malloc(strlen(xc_functional_keys[ii].name) + 1);
      strcpy(p,xc_functional_keys[ii].name);
      return p;
    }
  }
}


/*------------------------------------------------------*/
int xc_family_from_id(int id, int *family, int *number)
{
  int ii;

  /* first let us check if it is an LDA */
  for(ii=0; xc_lda_known_funct[ii]!=NULL; ii++){
    if(xc_lda_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_LDA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_LDA;
    }
  }

  /* or is it a hybrid LDA? */
  for(ii=0; xc_hyb_lda_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_lda_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_LDA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_LDA;
    }
  }

  /* or is it a GGA? */
  for(ii=0; xc_gga_known_funct[ii]!=NULL; ii++){
    if(xc_gga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_GGA;
    }
  }

  /* or is it a hybrid GGA? */
  for(ii=0; xc_hyb_gga_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_gga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_GGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_GGA;
    }
  }

  /* or is it a meta GGA? */
  for(ii=0; xc_mgga_known_funct[ii]!=NULL; ii++){
    if(xc_mgga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_MGGA;
    }
  }

  /* or is it a hybrid meta GGA? */
  for(ii=0; xc_hyb_mgga_known_funct[ii]!=NULL; ii++){
    if(xc_hyb_mgga_known_funct[ii]->number == id){
      if(family != NULL) *family = XC_FAMILY_HYB_MGGA;
      if(number != NULL) *number = ii;
      return XC_FAMILY_HYB_MGGA;
    }
  }

  return XC_FAMILY_UNKNOWN;
}

/*------------------------------------------------------*/
int xc_number_of_functionals()
{
  int num;

  for(num=0;;num++){
    if(xc_functional_keys[num].number == -1)
      return num;
  }

  fprintf(stderr, "Internal error in functionals.c\n");
  exit(1);
}

int xc_maximum_name_length()
{
  int i, N, maxlen, tmp;

  N=xc_number_of_functionals();

  maxlen=0;
  for(i=0;i<N;i++){
    tmp=strlen(xc_functional_keys[i].name);
    if(tmp > maxlen) maxlen=tmp;
  }

  return maxlen;
}

/*------------------------------------------------------*/
void xc_available_functional_numbers(int *list)
{
  int ii, N;
  N=xc_number_of_functionals();
  for(ii=0;ii<N;ii++){
    list[ii]=xc_functional_keys[ii].number;
  }
}

void xc_available_functional_names(char **list)
{
  int ii, N;

  N=xc_number_of_functionals();
  for(ii=0;ii<N;ii++) {
    strcpy(list[ii],xc_functional_keys[ii].name);
  }
}

/*------------------------------------------------------*/
xc_func_type *xc_func_alloc()
{
  xc_func_type *func;

  func = (xc_func_type *) libxc_malloc (sizeof (xc_func_type));
  return func;
}

/*------------------------------------------------------*/
int xc_func_init(xc_func_type *func, int functional, int nspin)
{
  int number;

  assert(func != NULL);
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);

  /* initialize structure */
  func->nspin       = nspin;

  func->params     = NULL;

  func->n_func_aux = 0;
  func->func_aux   = NULL;
  func->mix_coef   = NULL;
  func->cam_omega = func->cam_alpha = func->cam_beta = 0.0;
  func->nlc_b = func->nlc_C = 0.0;

  // we have to make a copy because the *_known_funct arrays live in
  // host memory (libxc_malloc instead returns memory than can be read
  // from GPU and CPU).
  xc_func_info_type * finfo = (xc_func_info_type *) libxc_malloc(sizeof(xc_func_info_type));
  
  switch(xc_family_from_id(functional, NULL, &number)){
  case(XC_FAMILY_LDA):
    *finfo = *xc_lda_known_funct[number];
    internal_counters_set_lda(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_LDA):
    *finfo = *xc_hyb_lda_known_funct[number];
    internal_counters_set_lda(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_GGA):
    *finfo = *xc_gga_known_funct[number];
    internal_counters_set_gga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_GGA):
    *finfo = *xc_hyb_gga_known_funct[number];
    internal_counters_set_gga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_MGGA):
    *finfo = *xc_mgga_known_funct[number];
    internal_counters_set_mgga(func->nspin, &(func->dim));
    break;

  case(XC_FAMILY_HYB_MGGA):
    *finfo = *xc_hyb_mgga_known_funct[number];
    internal_counters_set_mgga(func->nspin, &(func->dim));
    break;

  default:
    return -2; /* family not found */
  }

  func->info = finfo;
  
  /* see if we need to initialize the functional */
  if(func->info->init != NULL)
    func->info->init(func);

  /* see if we need to initialize the external parameters */
  if(func->info->ext_params.n > 0)
    func->info->ext_params.set(func, NULL);

  func->dens_threshold = func->info->dens_threshold;

  return 0;
}


/*------------------------------------------------------*/
void xc_func_end(xc_func_type *func)
{
  assert(func != NULL && func->info != NULL);

  /* call internal termination routine */
  if(func->info->end != NULL)
    func->info->end(func);

  /* terminate any auxiliary functional */
  if(func->n_func_aux > 0){
    int ii;

    for(ii=0; ii<func->n_func_aux; ii++){
      xc_func_end(func->func_aux[ii]);
      libxc_free(func->func_aux[ii]);
    }
    libxc_free(func->func_aux);
    func->n_func_aux = 0;
  }

  if(func->mix_coef != NULL){
    libxc_free(func->mix_coef);
    func->mix_coef = NULL;
  }

  /* deallocate any used parameter */
  if(func->params != NULL){
    libxc_free(func->params);
    func->params = NULL;
  }

  libxc_free((void *) func->info);
  func->info = NULL;  
}

/*------------------------------------------------------*/
void  xc_func_free(xc_func_type *p)
{
  libxc_free(p);
}

/*------------------------------------------------------*/
const xc_func_info_type *xc_func_get_info(const xc_func_type *p)
{
  return p->info;
}

/*------------------------------------------------------*/
void xc_func_set_dens_threshold(xc_func_type *p, double dens_threshold)
{
  int ii;

  p->dens_threshold = dens_threshold;

  for(ii=0; ii<p->n_func_aux; ii++) {
    xc_func_set_dens_threshold(p->func_aux[ii], dens_threshold);
  }
}

/*------------------------------------------------------*/
/* get/set external parameters                          */
void
xc_func_set_ext_params(xc_func_type *p, double *ext_params)
{
  assert(p->info->ext_params.n > 0);
  p->info->ext_params.set(p, ext_params);
}

void
xc_func_set_ext_params_name(xc_func_type *p, const char *name, double par)
{
  int ii;
  double *ext_params;

  assert(p != NULL && p->info->ext_params.n > 0);
  
  ext_params = libxc_malloc(p->info->ext_params.n*sizeof(double));
  for(ii=0; ii<p->info->ext_params.n; ii++){
    if(strcmp(p->info->ext_params.names[ii], name) == 0)
      ext_params[ii] = par;
    else
      ext_params[ii] = XC_EXT_PARAMS_DEFAULT;
  }
  xc_func_set_ext_params(p, ext_params);
  libxc_free(ext_params);
}


/*------------------------------------------------------*/
/* returns the mixing coefficient for the hybrid GGAs */
double xc_hyb_exx_coef(const xc_func_type *p)
{
  assert(p!=NULL);
 
  return p->cam_alpha;
}

/* returns the CAM parameters for screened hybrids */
void xc_hyb_cam_coef(const xc_func_type *p, double *omega, double *alpha, double *beta)
{
  assert(p!=NULL);

  *omega = p->cam_omega;
  *alpha = p->cam_alpha;
  *beta  = p->cam_beta;
}

/* returns the NLC parameters */
void xc_nlc_coef(const xc_func_type *p, double *nlc_b, double *nlc_C)
{
  assert(p!=NULL);

  *nlc_b = p->nlc_b;
  *nlc_C = p->nlc_C;
}
