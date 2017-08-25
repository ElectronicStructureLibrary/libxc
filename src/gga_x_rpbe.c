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

#define XC_GGA_X_RPBE  117 /* Hammer, Hansen & Norskov (PBE-like) */


typedef struct{
  double rpbe_kappa, rpbe_mu;
} gga_x_rpbe_params;


static void 
gga_x_rpbe_init(XC(func_type) *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = malloc(sizeof(gga_x_rpbe_params));

  /* same parameters as standard PBE */
  XC(gga_x_rpbe_set_params)(p, 0.8040, 0.2195149727645171);
}


void 
XC(gga_x_rpbe_set_params)(XC(func_type) *p, double kappa, double mu)
{
  gga_x_rpbe_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_x_rpbe_params *) (p->params);

  params->rpbe_kappa = kappa;
  params->rpbe_mu    = mu;
}

#include "maple2c/gga_x_rpbe.c"

#define func XC(gga_x_rpbe_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_rpbe) = {
  XC_GGA_X_RPBE,
  XC_EXCHANGE,
  "Hammer, Hansen, and Norskov",
  XC_FAMILY_GGA,
  {&xc_ref_Hammer1999_7413, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32,
  0, NULL, NULL,
  gga_x_rpbe_init, NULL, 
  NULL, work_gga_x, NULL
};
