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

#define XC_GGA_K_GE2          500 /* Second-order gradient expansion (integrated by parts) */
#define XC_GGA_K_VW           504 /* von Weiszaecker correction to Thomas-Fermi */

static void 
gga_k_ge2_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  /* value of beta in standard Becke 88 functional */
  switch(p->info->number){
  case XC_GGA_K_VW:
    p->func = 1; break;
  default: /* XC_GGA_K_GE2 */
    p->func = 0; break;
  }

static inline void
func(const XC(gga_type) *p, int order, FLOAT x,
     FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  FLOAT cnst[2] = {1.0/72.0, 1.0/8.0};

  *f = 1.0 + x*x/(cnst[p->func]*K_FACTOR_C);

  if(order < 1) return;

  *dfdx = 2.0*x/(cnst[p->func]*K_FACTOR_C);
  *ldfdx= 1.0  /(cnst[p->func]*K_FACTOR_C);
  
  if(order < 2) return;

  *d2fdx2 = 2.0/(cnst[p->func]*K_FACTOR_C);
}

#define XC_KINETIC_FUNCTIONAL
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_k_vw) = {
  XC_GGA_K_GE2,
  XC_KINETIC,
  "von Weiszaecker correction to Thomas-Fermi",
  XC_FAMILY_GGA,
  "CF von Weiszaecker, Z. Phys. 96, 431 (1935)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_k_ge2_init, 
  NULL, NULL,
  work_gga_x
};

const XC(func_info_type) XC(func_info_gga_k_ge2) = {
  XC_GGA_K_GE2,
  XC_KINETIC,
  "Second-order gradient expansion of the kinetic energy density",
  XC_FAMILY_GGA,
  "AS Kompaneets and ES Pavlovskii, Zh. Eksp. Teor. Fiz. 31, 427 (1956) [Sov. Phys. JETP 4, 328 (1957)]"
  "DA Kirznits, Zh. Eksp. Teor. Fiz. 32, 115 (1957) [Sov. Phys. JETP 5, 64 (1957)]",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  gga_k_ge2_init,
  NULL, NULL,
  work_gga_x
};
