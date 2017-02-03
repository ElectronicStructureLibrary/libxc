/*
 Copyright (C) 2014 Orbital-free DFT group at University of Florida, USA

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

/***********************************************************************
  Exchange and correlation free energy density and potential as 
  parametrized by 
    Valentin V. Karasiev, Travis Sjostrom, James Dufty, and S. B. Trickey
  Ported to C and libxc by Lazaro Calderin and Miguel Marques
************************************************************************/

#define XC_LDA_XC_KSDT    259    /* Karasiev et al. parametrization */

static const FLOAT ksdt_a[6] =
  {0.750, 3.043630, -0.0922700, 1.703500, 8.310510, 5.11050};
static const FLOAT ksdt_b[2][5] = {               /* b5 = Sqrt[3/2]/(lambda)*b3 */
  {0.2839970,  48.9321540, 0.3709190, 61.0953570, 0.871837422702767684673873513724},
  {0.3290010, 111.5983080, 0.5370530,105.0866630, 1.26233194679913807935662124247}
};
static const FLOAT ksdt_c[2][3] = {
  {0.8700890, 0.1930770, 2.4146440},
  {0.8489300, 0.1679520, 0.0888200}
};
static const FLOAT ksdt_d[2][5] = {
  {0.5798240,  94.5374540,  97.8396030,  59.9399990, 24.3880370},
  {0.5513300, 180.2131590, 134.4862310, 103.8616950, 17.7507100}
};
static const FLOAT ksdt_e[2][5] = {
  {0.2120360, 16.7312490, 28.4857920,  34.0288760, 17.2355150},
  {0.1531240, 19.5439450, 43.4003370, 120.2551450, 15.6628360}
};

typedef struct{
  FLOAT T;
} lda_xc_ksdt_params;


static void 
lda_xc_ksdt_init(XC(func_type) *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(lda_xc_ksdt_params));
}

#include "maple2c/lda_xc_ksdt.c"

#define func maple2c_func
#include "work_lda.c"

static const func_params_type ext_params[] = {
  {0.0, "Temperature"},
};

static void 
set_ext_params(XC(func_type) *p, const double *ext_params)
{
  lda_xc_ksdt_params *params;
  double ff;

  assert(p != NULL && p->params != NULL);
  params = (lda_xc_ksdt_params *) (p->params);

  ff = (ext_params == NULL) ? p->info->ext_params[0].value : ext_params[0];
  params->T = ff;
}

const XC(func_info_type) XC(func_info_lda_xc_ksdt) = {
  XC_LDA_XC_KSDT,
  XC_EXCHANGE_CORRELATION,
  "Karasiev, Sjostrom, Dufty & Trickey",
  XC_FAMILY_LDA,
  {&xc_ref_Karasiev2014_076403, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 0.0, 0.0, 1e-32,
  1, ext_params, set_ext_params,
  lda_xc_ksdt_init, NULL,
  work_lda, NULL, NULL
};
