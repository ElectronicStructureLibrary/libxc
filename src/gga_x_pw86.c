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
#include <assert.h>
#include "util.h"

#define XC_GGA_X_PW86         108 /* Perdew & Wang 86 */
#define XC_GGA_X_RPW86        144 /* refitted Perdew & Wang 86 */
#define XC_GGA_K_FR_PW86      515 /* Fuentealba & Reyes (PW86 version) */

typedef struct{
  FLOAT aa, bb, cc;
} gga_x_pw86_params;

static void 
gga_x_pw86_init(XC(func_type) *p)
{
  switch(p->info->number){
  case XC_GGA_X_RPW86:      p->func = 1; break;
  case XC_GGA_K_FR_PW86:    p->func = 2; break;
  default:                  p->func = 0; /* original PW86 */
  }
}

inline void
XC(gga_x_pw86_enhance)(const XC(func_type) *p, int order, FLOAT x, 
		       FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  static const gga_x_pw86_params par[3] = {
    {    1.296, 14.0,  0.2},
    {15*0.1234, 17.33, 0.163},
    {    2.208, 9.27,  0.2}
  };
  FLOAT ss, ss2, ss4, DD, dDD, d2DD, d3DD, aux;

  /* sanity check: do we have a valid functional */
  assert(p->func==0 || p->func==1 || p->func==2);

  ss     = X2S*x;
  ss2    = ss*ss;
  ss4    = ss2*ss2;

  DD     = 1.0 + par[p->func].aa*ss2 + par[p->func].bb*ss4 + par[p->func].cc*ss4*ss2;
  aux    = POW(DD, 1.0/15.0);
  *f     = aux;

  if(order < 1) return;

  dDD   = ss*(2.0*par[p->func].aa + 4.0*par[p->func].bb*ss2 + 6.0*par[p->func].cc*ss4);

  *dfdx  = X2S*dDD*aux/(15.0*DD);

  if(order < 2) return;

  d2DD    = 2.0*par[p->func].aa + 4.0*3.0*par[p->func].bb*ss2 + 6.0*5.0*par[p->func].cc*ss4;
  *d2fdx2 = X2S*X2S*aux*(15.0*DD*d2DD - 14.0*dDD*dDD)/(225.0*DD*DD);

  if(order < 3) return;

  d3DD    =  4.0*3.0*2.0*par[p->func].bb*ss + 6.0*5.0*4.0*par[p->func].cc*ss*ss2;
  *d3fdx3 = X2S*X2S*X2S*aux*(406.0*dDD*dDD*dDD - 630.0*DD*dDD*d2DD + 225.0*DD*DD*d3DD)/(3375.0*DD*DD*DD);
}

#define func XC(gga_x_pw86_enhance)
#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pw86) = {
  XC_GGA_X_PW86,
  XC_EXCHANGE,
  "Perdew & Wang 86",
  XC_FAMILY_GGA,
  "JP Perdew and Y Wang, Phys. Rev. B 33, 8800 (1986)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_pw86_init, NULL, NULL,
  work_gga_x,
  NULL
};

const XC(func_info_type) XC(func_info_gga_x_rpw86) = {
  XC_GGA_X_RPW86,
  XC_EXCHANGE,
  "Refitted Perdew & Wang 86",
  XC_FAMILY_GGA,
  "ED Murray, K Lee and DC Langreth, J. Chem. Theory Comput. 5, 2754-2762 (2009)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
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
  "P Fuentealba and O Reyes, Chem. Phys. Lett. 232, 31-34 (1995)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-32, 1e-32, 0.0, 1e-32,
  gga_x_pw86_init, NULL, NULL,
  work_gga_k,
  NULL
};
