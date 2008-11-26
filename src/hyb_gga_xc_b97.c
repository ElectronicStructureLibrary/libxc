/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_HYB_GGA_XC_B97      407 /* Becke 97                                 */
#define XC_HYB_GGA_XC_B97_1    408 /* Becke 97-1                               */
#define XC_HYB_GGA_XC_B97_2    410 /* Becke 97-2                               */
#define XC_HYB_GGA_XC_B97_K    413 /* Boese-Martin for Kinetics                */
#define XC_HYB_GGA_XC_B97_3    414 /* Boese-Martin for Kinetics                */

static void
hyb_gga_xc_b97_init(void *p_)
{
  const int iGGA[] = {XC_GGA_XC_B97, XC_GGA_XC_B97_1, XC_GGA_XC_B97_2, XC_GGA_XC_B97_K, XC_GGA_XC_B97_3};
  const FLOAT a0[] = {       0.1943,            0.21,            0.21,            0.42,    2.692880E-01};
  int func;

  XC(hyb_gga_type) *p = (XC(hyb_gga_type) *)p_;

  switch(p->info->number){
  case XC_HYB_GGA_XC_B97:      func = 0; break;
  case XC_HYB_GGA_XC_B97_1:    func = 1; break;
  case XC_HYB_GGA_XC_B97_2:    func = 2; break;
  case XC_HYB_GGA_XC_B97_K:    func = 3; break;
  case XC_HYB_GGA_XC_B97_3:    func = 4; break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_b97_init\n");
    exit(1);
    break;
  }

  p->lda_n = 0;
  p->gga_n = 1;

  XC(hyb_gga_alloc)(p);

  p->exx_coef = a0[func];

  XC(gga_init)(p->gga_aux[0], iGGA[func], p->nspin);
  p->gga_coef[0] = 1.0;
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_b1wc) = {
  XC_HYB_GGA_XC_B97,
  XC_EXCHANGE_CORRELATION,
  "Becke 97",
  XC_FAMILY_GGA,
  "AD Becke, J. Chem. Phys. 107, 8554 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  hyb_gga_xc_b97_init,
  NULL, NULL, NULL /* this is taken care by the generic routine */
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_1) = {
  XC_HYB_GGA_XC_B97_1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-1",
  XC_FAMILY_GGA,
  "FA Hamprecht, AJ Cohen, DJ Tozer, and NC Handy, J. Chem. Phys. 109, 6264 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  hyb_gga_xc_b97_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_2) = {
  XC_HYB_GGA_XC_B97_2,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-2",
  XC_FAMILY_GGA,
  "PJ Wilson, TJ Bradley, and DJ Tozer, J. Chem. Phys. 115, 9233 (2001)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_k) = {
  XC_HYB_GGA_XC_B97_K,
  XC_EXCHANGE_CORRELATION,
  "Boese-Martin for Kinetics",
  XC_FAMILY_GGA,
  "AD Boese and JML Martin, J. Chem. Phys., Vol. 121, 3405 (2004)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_3) = {
  XC_HYB_GGA_XC_B97_3,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-3",
  XC_FAMILY_GGA,
  "TW Keal and DJ Tozer, J. Chem. Phys. 123, 121103 (2005)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL
};
