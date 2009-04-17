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

#define XC_HYB_GGA_XC_B1WC   412 /* Becke 1-parameter mixture of WC and PBE */
#define XC_HYB_GGA_XC_B1LYP  416 /* Becke 1-parameter mixture of B88 and LYP */
#define XC_HYB_GGA_XC_B1PW91 417 /* Becke 1-parameter mixture of B88 and PW91 */
#define XC_HYB_GGA_XC_mPW1PW 418 /* Becke 1-parameter mixture of mPW91 and PW91 */
#define XC_HYB_GGA_XC_mPW1K  405 /* mixture of mPW91 and PW91 optimized for kinetics */

static void
hyb_gga_xc_b1_init(void *p_)
{
  const struct { FLOAT a0; int fx, fc; } 
  par[] = {
    {0.16,  XC_GGA_X_WC,    XC_GGA_C_PBE},  /* B1WC   */
    {0.25,  XC_GGA_X_B88,   XC_GGA_C_LYP},  /* B1LYP  */
    {0.25,  XC_GGA_X_B88,   XC_GGA_C_PW91}, /* B1PW91 */
    {0.25,  XC_GGA_X_mPW91, XC_GGA_C_PW91}, /* mPW1PW */
    {0.428, XC_GGA_X_mPW91, XC_GGA_C_PW91}  /* mPW1K  */
  };
  int func;

  XC(hyb_gga_type) *p = (XC(hyb_gga_type) *)p_;

  switch(p->info->number){
  case XC_HYB_GGA_XC_B1WC:   func = 0; break;
  case XC_HYB_GGA_XC_B1LYP:  func = 1; break;
  case XC_HYB_GGA_XC_B1PW91: func = 2; break;
  case XC_HYB_GGA_XC_mPW1PW: func = 3; break;
  case XC_HYB_GGA_XC_mPW1K:  func = 4; break;
  default:
    fprintf(stderr, "Internal error in gga_b97\n");
    exit(1);
    break;
  }

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  p->exx_coef = par[func].a0;

  XC(gga_init)(&p->mix->gga_mix[0], par[func].fx, p->nspin);
  p->mix->gga_coef[0] = (1.0 - par[func].a0);
  XC(gga_init)(&p->mix->gga_mix[1], par[func].fc, p->nspin);
  p->mix->gga_coef[1] = 1.0;
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_b1wc) = {
  XC_HYB_GGA_XC_B1WC,
  XC_EXCHANGE_CORRELATION,
  "B1WC",
  XC_FAMILY_HYB_GGA,
  "DI Bilc, R Orlando, R Shaltaf, G-M Rignanese, J Iniguez, and Ph Ghosez, Phys. Rev. B 77, 165107 (2008)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  hyb_gga_xc_b1_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b1lyp) = {
  XC_HYB_GGA_XC_B1LYP,
  XC_EXCHANGE_CORRELATION,
  "B1LYP",
  XC_FAMILY_HYB_GGA,
  "C. Adamo, V. Barone, Chem. Phys. Lett. 274, 242 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  hyb_gga_xc_b1_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b1pw91) = {
  XC_HYB_GGA_XC_B1PW91,
  XC_EXCHANGE_CORRELATION,
  "B1PW91",
  XC_FAMILY_HYB_GGA,
  "C. Adamo, V. Barone, Chem. Phys. Lett. 274, 242 (1997)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  hyb_gga_xc_b1_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw1pw) = {
  XC_HYB_GGA_XC_mPW1PW,
  XC_EXCHANGE_CORRELATION,
  "mPW1PW",
  XC_FAMILY_HYB_GGA,
  "C. Adamo, V. Barone, J. Chem. Phys. 108, 664 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  hyb_gga_xc_b1_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw1k) = {
  XC_HYB_GGA_XC_mPW1K,
  XC_EXCHANGE_CORRELATION,
  "mPW1K",
  XC_FAMILY_HYB_GGA,
  "BJ Lynch, PL Fast, M Harris, DGJ Truhlar, Phys. Chem. A 104, 4811 (2000)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  hyb_gga_xc_b1_init,
  NULL, NULL, NULL
};
