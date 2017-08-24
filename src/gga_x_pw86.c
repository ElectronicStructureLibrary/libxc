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

#define XC_GGA_X_PW86         108 /* Perdew & Wang 86 */
#define XC_GGA_X_RPW86        144 /* refitted Perdew & Wang 86 */
#define XC_GGA_K_FR_PW86      515 /* Fuentealba & Reyes (PW86 version) */

typedef struct{
  FLOAT aa, bb, cc;
} gga_x_pw86_params;


static const gga_x_pw86_params par_pw86 = {
  1.296, 14.0, 0.2
};

static const gga_x_pw86_params par_rpw86 = {
  15.0*0.1234, 17.33, 0.163,
};

static const gga_x_pw86_params par_fr_pw86 = {
  2.208, 9.27, 0.2
};

static void 
gga_x_pw86_init(XC(func_type) *p)
{
  gga_x_pw86_params *params;

  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_pw86_params));
  params = (gga_x_pw86_params *) (p->params);

  switch(p->info->number){
  case XC_GGA_X_PW86: 
    memcpy(params, &par_pw86, sizeof(gga_x_pw86_params));
    break;
  case XC_GGA_X_RPW86:
    memcpy(params, &par_rpw86, sizeof(gga_x_pw86_params));
    break;
  case XC_GGA_K_FR_PW86:
    memcpy(params, &par_fr_pw86, sizeof(gga_x_pw86_params));
    break;
  default:
    fprintf(stderr, "Internal error in gga_x_pw86\n");
    exit(1);
  }
}

#include "maple2c/gga_x_pw86.c"

#define func maple2c_func
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pw86) = {
  XC_GGA_X_PW86,
  XC_EXCHANGE,
  "Perdew & Wang 86",
  XC_FAMILY_GGA,
  {&xc_ref_Perdew1986_8800, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0,
  0, NULL, NULL,
  gga_x_pw86_init, NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_rpw86) = {
  XC_GGA_X_RPW86,
  XC_EXCHANGE,
  "Refitted Perdew & Wang 86",
  XC_FAMILY_GGA,
  {&xc_ref_Murray2009_2754, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0,
  0, NULL, NULL,
  gga_x_pw86_init, NULL, NULL,
  work_gga_x,
  NULL
};

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_fr_pw86) = {
  XC_GGA_K_FR_PW86,
  XC_KINETIC,
  "Fuentealba & Reyes (PW86 version)",
  XC_FAMILY_GGA,
  {&xc_ref_Fuentealba1995_31, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0,
  0, NULL, NULL,
  gga_x_pw86_init, NULL, NULL,
  work_gga_k,
  NULL
};
