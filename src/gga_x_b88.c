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

#include "util.h"

#define XC_GGA_X_B88          106 /* Becke 88 */
#define XC_GGA_X_OPTB88_VDW   139 /* Becke 88 reoptimized to be used with vdW functional of Dion et al */
#define XC_GGA_X_MB88         149 /* Modified Becke 88 for proton transfer */
#define XC_GGA_X_EB88         271 /* Non-empirical (excogitated) B88 functional of Becke and Elliott */
#define XC_GGA_K_LLP          522 /* Lee, Lee & Parr */
#define XC_GGA_K_FR_B88       514 /* Fuentealba & Reyes (B88 version) */
#define XC_GGA_X_B88M         570 /* Becke 88 reoptimized to be used with mgga_c_tau1 */

typedef struct{
  double beta, gamma;
} gga_x_b88_params;


static void 
gga_x_b88_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_b88_params));

  /* value of beta in standard Becke 88 functional */
  switch(p->info->number){
  case XC_GGA_X_B88:
    xc_gga_x_b88_set_params(p, 0.0042, 6.0);
    break;
  case XC_GGA_X_OPTB88_VDW:
    xc_gga_x_b88_set_params(p, 0.00336865923905927, 6.98131700797731);
    break;
  case XC_GGA_K_LLP:
    xc_gga_x_b88_set_params(p, X_FACTOR_C*0.0044188, 0.0253/(X_FACTOR_C*0.0044188));
    break;
  case XC_GGA_K_FR_B88:
    xc_gga_x_b88_set_params(p, X_FACTOR_C*0.004596, 0.02774/(X_FACTOR_C*0.004596));
    break;
  case XC_GGA_X_MB88:
    xc_gga_x_b88_set_params(p, 0.0011, 6.0);
    break;
  case XC_GGA_X_EB88:
    xc_gga_x_b88_set_params(p, 0.0050/M_CBRT2, 6.0);
    break;
  case XC_GGA_X_B88M:
    xc_gga_x_b88_set_params(p, 0.0045, 6.0);
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_b88\n");
    exit(1);
  }
}


void 
xc_gga_x_b88_set_params(xc_func_type *p, double beta, double gamma)
{
  gga_x_b88_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_b88_params *) (p->params);

  params->beta  = beta;
  params->gamma = gamma;
}


#include "maple2c/gga_x_b88.c"

#define func xc_gga_x_b88_enhance
#include "work_gga_x.c"

const xc_func_info_type xc_func_info_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1988_3098, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL, 
  NULL, work_gga_x, NULL
};

const xc_func_info_type xc_func_info_gga_x_optb88_vdw = {
  XC_GGA_X_OPTB88_VDW,
  XC_EXCHANGE,
  "opt-Becke 88 for vdW",
  XC_FAMILY_GGA,
  {&xc_ref_Klimes2010_022201, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL, 
  NULL, work_gga_x, NULL
};

const xc_func_info_type xc_func_info_gga_x_mb88 = {
  XC_GGA_X_MB88,
  XC_EXCHANGE,
  "Modified Becke 88 for proton transfer",
  XC_FAMILY_GGA,
  {&xc_ref_Tognetti2009_14415, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL, 
  NULL, work_gga_x, NULL
};

const xc_func_info_type xc_func_info_gga_x_eb88 = {
  XC_GGA_X_EB88,
  XC_EXCHANGE,
  "Non-empirical (excogitated) B88 functional of Becke and Elliott",
  XC_FAMILY_GGA,
  {&xc_ref_Elliott2009_1485, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init,  NULL, 
  NULL, work_gga_x, NULL
};

const xc_func_info_type xc_func_info_gga_x_b88m = {
  XC_GGA_X_B88M,
  XC_EXCHANGE,
  "Becke 88 reoptimized to be used with tau1",
  XC_FAMILY_GGA,
  {&xc_ref_Proynov2000_10013, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init,  NULL, 
  NULL, work_gga_x, NULL
};

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const xc_func_info_type xc_func_info_gga_k_llp = {
  XC_GGA_K_LLP,
  XC_KINETIC,
  "Lee, Lee & Parr",
  XC_FAMILY_GGA,
  {&xc_ref_Lee1991_768, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL,
  NULL, work_gga_k, NULL
};

const xc_func_info_type xc_func_info_gga_k_fr_b88 = {
  XC_GGA_K_FR_B88,
  XC_KINETIC,
  "Fuentealba & Reyes (B88 version)",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-25,
  0, NULL, NULL,
  gga_x_b88_init, NULL,
  NULL, work_gga_k, NULL
};
