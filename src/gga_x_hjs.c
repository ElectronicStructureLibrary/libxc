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

#define XC_GGA_X_HJS_PBE     525 /* HJS screened exchange PBE version */
#define XC_GGA_X_HJS_PBE_SOL 526 /* HJS screened exchange PBE_SOL version */
#define XC_GGA_X_HJS_B88     527 /* HJS screened exchange B88 version */
#define XC_GGA_X_HJS_B97X    528 /* HJS screened exchange B97x version */

typedef struct{
  FLOAT omega;

  const FLOAT *a, *b; /* pointers to the a and b parameters */
} gga_x_hjs_params;

static const FLOAT a_PBE[] = 
  {0.0159941, 0.0852995, -0.160368, 0.152645, -0.0971263, 0.0422061};
static const FLOAT b_PBE[] = 
  {5.33319, -12.4780, 11.0988, -5.11013, 1.71468, -0.610380, 0.307555, -0.0770547, 0.0334840};

static const FLOAT a_PBE_sol[] = 
  {0.0047333, 0.0403304, -0.0574615, 0.0435395, -0.0216251, 0.0063721};
static const FLOAT b_PBE_sol[] = 
  {8.52056, -13.9885, 9.28583, -3.27287, 0.843499, -0.235543, 0.0847074, -0.0171561, 0.0050552};

static const FLOAT a_B88[] =
  {0.00968615, -0.0242498, 0.0259009, -0.0136606, 0.00309606, -7.32583e-5};
static const FLOAT b_B88[] =
  {-2.50356, 2.79656, -1.79401, 0.714888, -0.165924, 0.0118379, 0.0037806, -1.57905e-4, 1.45323e-6};

static const FLOAT a_B97x[] =
  {0.0027355, 0.0432970, -0.0669379, 0.0699060, -0.0474635, 0.0153092};
static const FLOAT b_B97x[] =
  {15.8279, -26.8145, 17.8127, -5.98246, 1.25408, -0.270783, 0.0919536, -0.0140960, 0.0045466};

static const FLOAT a_B88_V2[] =
  {0.0253933, -0.0673075, 0.0891476, -0.0454168, -0.00765813, 0.0142506};
static const FLOAT b_B88_V2[] =
  {-2.6506, 3.91108, -3.31509, 1.54485, -0.198386, -0.136112, 0.0647862, 0.0159586, -0.000245066};

static void
gga_x_hjs_init(XC(func_type) *p)
{
  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_x_hjs_params));

  /* we take 0.11 as the default for hjs */
  XC(gga_x_hjs_set_params)(p, 0.11);

  switch(p->info->number){
  case XC_GGA_X_HJS_PBE:
    ((gga_x_hjs_params *)(p->params))->a = a_PBE;
    ((gga_x_hjs_params *)(p->params))->b = b_PBE;
    break;
  case XC_GGA_X_HJS_PBE_SOL:
    ((gga_x_hjs_params *)(p->params))->a = a_PBE_sol;
    ((gga_x_hjs_params *)(p->params))->b = b_PBE_sol;
    break;
  case XC_GGA_X_HJS_B88:
    ((gga_x_hjs_params *)(p->params))->a = a_B88;
    ((gga_x_hjs_params *)(p->params))->b = b_B88;
    break;
  case XC_GGA_X_HJS_B97X:
    ((gga_x_hjs_params *)(p->params))->a = a_B97x;
    ((gga_x_hjs_params *)(p->params))->b = b_B97x;
    break;
  case XC_GGA_X_HJS_B88_V2:
    ((gga_x_hjs_params *)(p->params))->a = a_B88_V2;
    ((gga_x_hjs_params *)(p->params))->b = b_B88_V2;
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_hjs_init\n");
    exit(1);
  }
}

void 
XC(gga_x_hjs_set_params)(XC(func_type) *p, FLOAT omega)
{
  gga_x_hjs_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_hjs_params *) (p->params);

  params->omega = omega;
}


#include "maple2c/gga_x_hjs.c"

#define func maple2c_func
#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_x_hjs_pbe) = {
  XC_GGA_X_HJS_PBE,
  XC_EXCHANGE,
  "HJS screened exchange PBE version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_hjs_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_x_hjs_pbe_sol) = {
  XC_GGA_X_HJS_PBE_SOL,
  XC_EXCHANGE,
  "HJS screened exchange PBE_SOL version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_hjs_init, NULL, 
  NULL, work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_x_hjs_b88) = {
  XC_GGA_X_HJS_B88,
  XC_EXCHANGE,
  "HJS screened exchange B88 version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-6, 1e-11, 0.0, 0.0, /* densities smaller than 1e-6 yield NaNs */
  0, NULL, NULL,
  gga_x_hjs_init, NULL, 
  NULL,  work_gga_c, NULL
};

const XC(func_info_type) XC(func_info_gga_x_hjs_b97x) = {
  XC_GGA_X_HJS_B97X,
  XC_EXCHANGE,
  "HJS screened exchange B97x version",
  XC_FAMILY_GGA,
  {&xc_ref_Henderson2008_194105, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  0, NULL, NULL,
  gga_x_hjs_init, NULL, 
  NULL, work_gga_c, NULL
};
