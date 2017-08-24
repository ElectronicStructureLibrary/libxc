/*
 Copyright (C) 2008 M.A.L. Marques

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

#define XC_MGGA_C_M06_L         233 /* M06-Local functional from Minnesota */
#define XC_MGGA_C_M06_HF        234 /* Worker for M06-HF functional        */
#define XC_MGGA_C_M06           235 /* Worker for M06 functional           */
#define XC_MGGA_C_M06_2X        236 /* Worker for M06-2X functional        */

typedef struct{
  FLOAT gamma_ss, gamma_ab, alpha_ss, alpha_ab;
  const FLOAT css[5], cab[5], dss[6], dab[6];
} mgga_c_m06l_params;

static const mgga_c_m06l_params par_m06l = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  { 5.349466e-01,  5.396620e-01, -3.161217e+01,  5.149592e+01, -2.919613e+01},
  { 6.042374e-01,  1.776783e+02, -2.513252e+02,  7.635173e+01, -1.255699e+01},
  { 4.650534e-01,  1.617589e-01,  1.833657e-01,  4.692100e-04, -4.990573e-03,  0.000000e+00},
  { 3.957626e-01, -5.614546e-01,  1.403963e-02,  9.831442e-04, -3.577176e-03,  0.000000e+00}
};

static const mgga_c_m06l_params par_m06hf = {
  0.06, 0.0031, 0.00515088, 0.00304966,
  { 1.023254e-01, -2.453783e+00,  2.913180e+01, -3.494358e+01,  2.315955e+01},
  { 1.674634e+00,  5.732017e+01,  5.955416e+01, -2.311007e+02,  1.255199e+02},
  { 8.976746e-01, -2.345830e-01,  2.368173e-01, -9.913890e-04, -1.146165e-02,  0.000000e+00},
  {-6.746338e-01, -1.534002e-01, -9.021521e-02, -1.292037e-03, -2.352983e-04,  0.000000e+00}
};

static const mgga_c_m06l_params par_m06 = {
  0.06, 0.0031, 0.00515088, 0.00304966,  
  { 5.094055e-01, -1.491085e+00,  1.723922e+01, -3.859018e+01,  2.845044e+01},
  { 3.741539e+00,  2.187098e+02, -4.531252e+02,  2.936479e+02, -6.287470e+01},
  { 4.905945e-01, -1.437348e-01,  2.357824e-01,  1.871015e-03, -3.788963e-03,  0.000000e+00},
  {-2.741539e+00, -6.720113e-01, -7.932688e-02,  1.918681e-03, -2.032902e-03,  0.000000e+00}
};

static const mgga_c_m06l_params par_m062x = {
  0.06, 0.0031, 0.00515088, 0.00304966,  
  { 3.097855e-01, -5.528642e+00,  1.347420e+01, -3.213623e+01,  2.846742e+01},
  { 8.833596e-01,  3.357972e+01, -7.043548e+01,  4.978271e+01, -1.852891e+01},
  { 6.902145e-01,  9.847204e-02,  2.214797e-01, -1.968264e-03, -6.775479e-03,  0.000000e+00},
  { 1.166404e-01, -9.120847e-02, -6.726189e-02,  6.720580e-05,  8.448011e-04,  0.000000e+00}
};

static void 
mgga_c_vsxc_init(XC(func_type) *p)
{
  mgga_c_m06l_params *params;

  assert(p != NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_C_PW_MOD, XC_POLARIZED);

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(mgga_c_m06l_params));
  params = (mgga_c_m06l_params *)p->params;

  switch(p->info->number){
  case XC_MGGA_C_M06_L:
    memcpy(params, &par_m06l, sizeof(mgga_c_m06l_params));
    break;
  case XC_MGGA_C_M06_HF: 
    memcpy(params, &par_m06hf, sizeof(mgga_c_m06l_params));
    break;
  case XC_MGGA_C_M06:
    memcpy(params, &par_m06, sizeof(mgga_c_m06l_params));
    break;
  case XC_MGGA_C_M06_2X:
    memcpy(params, &par_m062x, sizeof(mgga_c_m06l_params));
    break;
  default:
    fprintf(stderr, "Internal error in mgga_c_m06l\n");
    exit(1);
  }  
}

#include "maple2c/mgga_c_m06l.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_m06_l) = {
  XC_MGGA_C_M06_L,
  XC_CORRELATION,
  "Minnesota M06-L functional",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_194101, &xc_ref_Zhao2008_215, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU,
  0, NULL, NULL,
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga_c,
};

const XC(func_info_type) XC(func_info_mgga_c_m06_hf) = {
  XC_MGGA_C_M06_HF,
  XC_CORRELATION,
  "Worker for hyb_mgga_xc_m06_hf",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2006_13126, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU,
  0, NULL, NULL,
  mgga_c_vsxc_init, NULL, 
  NULL, NULL, work_mgga_c,
};

const XC(func_info_type) XC(func_info_mgga_c_m06) = {
  XC_MGGA_C_M06,
  XC_CORRELATION,
  "Worker for hyb_mgga_xc_m06",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU,
  0, NULL, NULL,
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga_c,
};

const XC(func_info_type) XC(func_info_mgga_c_m06_2x) = {
  XC_MGGA_C_M06_2X,
  XC_CORRELATION,
  "Worker for hyb_mgga_xc_m06_2x",
  XC_FAMILY_MGGA,
  {&xc_ref_Zhao2008_215, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  MIN_DENS, MIN_GRAD, MIN_TAU,
  0, NULL, NULL,
  mgga_c_vsxc_init, NULL,
  NULL, NULL, work_mgga_c,
};
