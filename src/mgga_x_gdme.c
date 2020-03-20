/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

/* Local tau approximation */

#define XC_MGGA_X_GDME_NV   687 /* Generalized density-matrix with a=1/2      */
#define XC_MGGA_X_GDME_0    689 /* Generalized density-matrix with a=0        */
#define XC_MGGA_X_GDME_KOS  690 /* Generalized density-matrix with a=0.00638  */
#define XC_MGGA_X_GDME_VT   691 /*  Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria */

typedef struct{
  double a, AA, BB;
} mgga_x_gdme_params;


static void 
mgga_x_gdme_init(xc_func_type *p)
{
  mgga_x_gdme_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_gdme_params));
  params = (mgga_x_gdme_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_GDME_NV:
    /* default values set by set_ext_params */
    break;
  case XC_MGGA_X_GDME_0:
    params->a  = 0.0;
    params->AA = 9.0*M_PI/4.0;
    params->BB = 35.0*M_PI/12.0;
    break;
  case XC_MGGA_X_GDME_KOS:
    params->a  = 0.00638;
    params->AA = 9.0*M_PI/4.0;
    params->BB = 35.0*M_PI/12.0;
    break;
  case XC_MGGA_X_GDME_VT:
    params->a  = 0.0;
    params->AA = 7.31275;
    params->BB = 5.43182;
    break;    
  default:
    fprintf(stderr, "Internal error in mgga_x_gdme\n");
    exit(1);
  }
}

static func_params_type ext_params[] = {
  {"_a",  0.5,            "center of the s expansion of density-matrix"},
  {"_AA", 9.0*M_PI/4.0,   "parameter of the first (LDA) term"},
  {"_BB", 35.0*M_PI/12.0, "parameter of the correction term"}
};

static void 
set_ext_params(xc_func_type *p, const double *ext_params)
{
  mgga_x_gdme_params *params;

  assert(p != NULL && p->params != NULL);
  params = (mgga_x_gdme_params *) (p->params);

  params->a  = get_ext_param(p->info->ext_params, ext_params, 0);
  params->AA = get_ext_param(p->info->ext_params, ext_params, 1);
  params->BB = get_ext_param(p->info->ext_params, ext_params, 2);
}



#include "decl_mgga.h"
#include "maple2c/mgga_exc/mgga_x_gdme.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_nv = {
  XC_MGGA_X_GDME_NV,
  XC_EXCHANGE,
  "Generalized density-matrix with a=1/2",
  XC_FAMILY_MGGA,
  {&xc_ref_Negele1972_1472, &xc_ref_Koehl1996_835, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  3, ext_params, set_ext_params,
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_0 = {
  XC_MGGA_X_GDME_0,
  XC_EXCHANGE,
  "Generalized density-matrix with a=0",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_kos = {
  XC_MGGA_X_GDME_KOS,
  XC_EXCHANGE,
  "Generalized density-matrix with a=0.00638",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_gdme_vt = {
  XC_MGGA_X_GDME_VT,
  XC_EXCHANGE,
  "Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria",
  XC_FAMILY_MGGA,
  {&xc_ref_Koehl1996_835, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | MAPLE2C_FLAGS,
  1.0e-23,
  {0, NULL, NULL, NULL, NULL},
  mgga_x_gdme_init, NULL,
  NULL, NULL, work_mgga,
};

