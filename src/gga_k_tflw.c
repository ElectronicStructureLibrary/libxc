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

/* for a review on the values of lambda and gamma, please see EV
Ludena and VV Karasiev, in "Reviews of Modern Quantum Chemistry: a
Celebration of the Contributions of Robert G. Parr, edited by KD Sen
(World Scientific, Singapore, 2002), p. 612.
 */

#define XC_GGA_K_TFVW          52  /* Thomas-Fermi plus von Weiszaecker correction */
#define XC_GGA_K_VW            500 /* von Weiszaecker functional */
#define XC_GGA_K_GE2           501 /* Second-order gradient expansion (l = 1/9) */
#define XC_GGA_K_GOLDEN        502 /* TF-lambda-vW form by Golden (l = 13/45) */
#define XC_GGA_K_YT65          503 /* TF-lambda-vW form by Yonei and Tomishima (l = 1/5) */
#define XC_GGA_K_BALTIN        504 /* TF-lambda-vW form by Baltin (l = 5/9) */
#define XC_GGA_K_LIEB          505 /* TF-lambda-vW form by Lieb (l = 0.185909191) */
#define XC_GGA_K_ABSP1         506 /* gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)] */
#define XC_GGA_K_ABSP2         507 /* gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)] */
#define XC_GGA_K_ABSP3         277 /* gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35] */
#define XC_GGA_K_ABSP4         278 /* gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))] */
#define XC_GGA_K_GR            508 /* gamma-TFvW form by Gazquez and Robles */
#define XC_GGA_K_LUDENA        509 /* gamma-TFvW form by Ludena */
#define XC_GGA_K_GP85          510 /* gamma-TFvW form by Ghosh and Parr */

typedef struct{
  FLOAT gamma, lambda;
} gga_k_tflw_params;


static void 
gga_k_tflw_init(XC(func_type) *p)
{

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_k_tflw_params));

  /* This automatically sets gamma and lambda depending on the functional chosen.
     We put by default N = 1.0 */
  XC(gga_k_tflw_set_params)(p, -1.0, -1.0, 1.0);
}

/* for automatically assigning lambda and gamma set them to -1 */
void 
XC(gga_k_tflw_set_params)(XC(func_type) *p, FLOAT gamma, FLOAT lambda, FLOAT N)
{
  gga_k_tflw_params *params;
  FLOAT C0 = CBRT(M_PI/3.0);
  FLOAT C1 = CBRT(M_PI*M_PI/36.0)/6.0 - CBRT(M_PI*M_PI/9.0)/4.0;
  
  assert(p != NULL && p->params != NULL);
  params = (gga_k_tflw_params *) (p->params);

  params->gamma = 1.0;
  if(gamma > 0.0){
    params->gamma = gamma;
  }else if(N > 0.0){
    switch(p->info->number){
    case XC_GGA_K_TFVW:
      params->gamma = 1.0;
      break;
    case XC_GGA_K_VW:
      params->gamma = 0.0;
      break;
    case XC_GGA_K_ABSP1:      /* Ref. 79 */
      params->gamma = 1.0 - 1.412/CBRT(N);
      break;
    case XC_GGA_K_ABSP2:      /* Ref. 79 */
      params->gamma = 1.0 - 1.332/CBRT(N);
      break;
    case XC_GGA_K_ABSP3:      /* Ref. 79 */
      params->gamma = 1.0 - 1.513/POW(N, 0.35);
      break;
    case XC_GGA_K_ABSP4:      /* Ref. 79 */
      params->gamma = 1.0/(1.0 + 1.332/CBRT(N));
      break;
    case XC_GGA_K_GR:         /* Ref. 80 */
      params->gamma = (1.0 - 2.0/N)*(1.0 - C0/CBRT(N) + C1*CBRT(N*N));
      break;
    case XC_GGA_K_LUDENA:     /* Ref. 82 */
      params->gamma = CBRT(6.0*M_PI)*M_PI*M_PI*(1.0 - 1.0/(N*N));
	break;
    case XC_GGA_K_GP85:       /* Ref. 86 */
      params->gamma = CBRT(6.0*M_PI*M_PI)*M_PI*M_PI/4.0*
	(1.0 - 1.0/N)*(1.0 + 1.0/N + 6.0/(N*N));
      break;
    }
  }

  params->lambda = 1.0;
  if(lambda > 0.0){
    params->lambda  = lambda;
  }else{
    switch(p->info->number){
    case XC_GGA_K_TFVW:
      params->lambda = 1.0;
      break;
    case XC_GGA_K_GE2:
      params->lambda = 1.0/9.0;
      break;
    case XC_GGA_K_GOLDEN:     /* Ref. 33 */
      params->lambda = 13.0/45.0;
      break;
    case XC_GGA_K_YT65:       /* Ref. 57 */
      params->lambda = 1.0/5.0;
      break;
    case XC_GGA_K_BALTIN:     /* Ref. 66 */
      params->lambda = 5.0/9.0;
      break;
    case XC_GGA_K_LIEB:       /* Ref. 12 */
      params->lambda = 0.185909191;   /* 1/5.37897... */
      break;
    case XC_GGA_K_ABSP4:      /* Ref. 79 */
      params->lambda = 1.0/(1.0 + 1.332/CBRT(N));
      break;
    }
  }
}


static inline void 
func(const XC(func_type) *p, int order, FLOAT x, 
     FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT lambda, gamma;

  assert(p->params != NULL);
  lambda  = ((gga_k_tflw_params *) (p->params))->lambda;
  gamma   = ((gga_k_tflw_params *) (p->params))->gamma;

  lambda /= 8.0; /* the von Weiszaecker coefficient */
  
  *f = gamma + lambda*x*x/K_FACTOR_C;

  if(order < 1) return;

  *dfdx = 2.0*lambda*x/K_FACTOR_C;
  
  if(order < 2) return;

  *d2fdx2 = 2.0*lambda/K_FACTOR_C;
}

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_tfvw) = {
  XC_GGA_K_TFVW,
  XC_KINETIC,
  "Thomas-Fermi plus von Weiszaecker correction",
  XC_FAMILY_GGA,
  {&xc_ref_Weizsacker1935_431, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init, 
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_vw) = {
  XC_GGA_K_VW,
  XC_KINETIC,
  "von Weiszaecker correction to Thomas-Fermi",
  XC_FAMILY_GGA,
  {&xc_ref_Weizsacker1935_431, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init, 
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_ge2) = {
  XC_GGA_K_GE2,
  XC_KINETIC,
  "Second-order gradient expansion of the kinetic energy density",
  XC_FAMILY_GGA,
  {&xc_ref_Kompaneets1956_427, &xc_ref_Kirznits1957_115, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_golden) = {
  XC_GGA_K_GOLDEN,
  XC_KINETIC,
  "TF-lambda-vW form by Golden (l = 13/45)",
  XC_FAMILY_GGA,
  {&xc_ref_Golden1957_604, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_yt65) = {
  XC_GGA_K_YT65,
  XC_KINETIC,
  "TF-lambda-vW form by Yonei and Tomishima (l = 1/5)",
  XC_FAMILY_GGA,
  {&xc_ref_Yonei1965_1051, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_baltin) = {
  XC_GGA_K_BALTIN,
  XC_KINETIC,
  "TF-lambda-vW form by Baltin (l = 5/9)",
  XC_FAMILY_GGA,
  {&xc_ref_Baltin1972_1176, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_lieb) = {
  XC_GGA_K_LIEB,
  XC_KINETIC,
  "TF-lambda-vW form by Lieb (l = 0.185909191)",
  XC_FAMILY_GGA,
  {&xc_ref_Lieb1981_603, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_absp1) = {
  XC_GGA_K_ABSP1,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_absp2) = {
  XC_GGA_K_ABSP2,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_absp3) = {
  XC_GGA_K_ABSP3,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_absp4) = {
  XC_GGA_K_ABSP4,
  XC_KINETIC,
  "gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))]",
  XC_FAMILY_GGA,
  {&xc_ref_Acharya1980_6978, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_gr) = {
  XC_GGA_K_GR,
  XC_KINETIC,
  "gamma-TFvW form by Gazquez and Robles",
  XC_FAMILY_GGA,
  {&xc_ref_Gazquez1982_1467, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_ludena) = {
  XC_GGA_K_LUDENA,
  XC_KINETIC,
  "gamma-TFvW form by Ludena",
  XC_FAMILY_GGA,
  {&xc_ref_Ludena1986, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};

const XC(func_info_type) XC(func_info_gga_k_gp85) = {
  XC_GGA_K_GP85,
  XC_KINETIC,
  "gamma-TFvW form by Ghosh and Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Ghosh1985_3307, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_k_tflw_init,
  NULL, NULL,
  work_gga_k,
  NULL
};
