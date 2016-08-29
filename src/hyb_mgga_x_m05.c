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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_HYB_MGGA_X_M05      438 /* M05 functional from Minnesota     */
#define XC_HYB_MGGA_X_M05_2X   439 /* M05-2X functional from Minnesota  */
#define XC_HYB_MGGA_X_M06_2X   450 /* M06-2X functional from Minnesota  */
#define XC_HYB_MGGA_X_DLDF      36 /* Dispersionless Density Functional */


static const FLOAT a_m05[12] = 
  {1.0, 0.08151, -0.43956, -3.22422, 2.01819, 8.79431, -0.00295,
   9.82029, -4.82351, -48.17574, 3.64802, 34.02248};

static const FLOAT a_m05_2x[12] =
  {1.0, -0.56833, -1.30057, 5.50070, 9.06402, -32.21075, -23.73298,
   70.22996, 29.88614, -60.25778, -13.22205, 15.23694};

static const FLOAT a_m06_2x[12] =
  {4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00, -1.422133e+01,
   1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,  1.522808e+01,  1.522227e+01};

static const FLOAT a_dldf[5] =
  {1.0, -0.1637571, -0.1880028, -0.4490609, -0.0082359};

typedef struct{
  int n;
  const FLOAT *a;
  FLOAT csi_HF;
} mgga_x_m05_params;


static void
mgga_x_m05_init(XC(func_type) *p)
{
  mgga_x_m05_params *params;

  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_X_PBE, p->nspin);

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m05_params));
  params = (mgga_x_m05_params *) (p->params);

  switch(p->info->number){
  case XC_HYB_MGGA_X_M05: 
    params->n      = 12;
    params->a      = a_m05;
    p->cam_alpha   = 0.28;
    params->csi_HF = 1.0 - p->cam_alpha;
    break;
  case XC_HYB_MGGA_X_M05_2X:
    params->n      = 12;
    params->a      = a_m05_2x;
    p->cam_alpha   = 0.56;
    params->csi_HF = 1.0 - p->cam_alpha;
    break;
  case XC_HYB_MGGA_X_M06_2X:
    params->n      = 12;
    params->a      = a_m06_2x;
    p->cam_alpha   = 0.54;
    params->csi_HF = 1.0; /* the mixing is already included in the params->a */
    break;
  case XC_HYB_MGGA_X_DLDF:
    params->n      = 5;
    params->a      = a_dldf;
    p->cam_alpha   = 0.6144129;
    params->csi_HF = 1.0 - p->cam_alpha;

    XC(gga_x_pbe_set_params)(p->func_aux[0], 4.8827323, 0.3511128);

    break;
  default:
    fprintf(stderr, "Internal error in mgga_x_m05\n");
    exit(1);
  }

}


static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  mgga_x_m05_params *params;

  FLOAT e_f, e_dfdx, e_d2fdx2;
  FLOAT fw, dfwdt;

  assert(pt != NULL && pt->params != NULL);
  params = (mgga_x_m05_params *) (pt->params);
  
  XC(gga_x_pbe_enhance)(pt->func_aux[0], r->order, r->x, &e_f, &e_dfdx, &e_d2fdx2, NULL);
  
  XC(mgga_series_w)(r->order, params->n, params->a, r->t, &fw, &dfwdt);

  r->f = params->csi_HF*e_f*fw;

  if(r->order < 1) return;

  r->dfdx = params->csi_HF*e_dfdx*fw;
  r->dfdt = params->csi_HF*e_f*dfwdt;

  if(r->order < 2) return;

}


#include "work_mgga_x.c"


const XC(func_info_type) XC(func_info_hyb_mgga_x_m05) = {
  XC_HYB_MGGA_X_M05,
  XC_EXCHANGE,
  "Minnesota M05 functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_161103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  mgga_x_m05_init,
  NULL, NULL, NULL,
  work_mgga_x,
};


const XC(func_info_type) XC(func_info_hyb_mgga_x_m05_2x) = {
  XC_HYB_MGGA_X_M05_2X,
  XC_EXCHANGE,
  "M05-2X functional from Minnesota",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2006_364, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  mgga_x_m05_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_hyb_mgga_x_m06_2x) = {
  XC_HYB_MGGA_X_M06_2X,
  XC_EXCHANGE,
  "M06-2X functional from Minnesota",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  0, NULL, NULL,
  mgga_x_m05_init,
  NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};

const XC(func_info_type) XC(func_info_hyb_mgga_x_dldf) = {
  XC_HYB_MGGA_X_DLDF,
  XC_EXCHANGE,
  "Dispersionless Density Functional",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Pernal2009_263201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU, MIN_ZETA,
  0, NULL, NULL,
  mgga_x_m05_init,
  NULL,
  NULL, NULL,        /* this is not an LDA                   */
  work_mgga_x,
};
