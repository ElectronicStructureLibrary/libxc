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

#define XC_HYB_GGA_XC_B3PW91  401 /* The original hybrid proposed by Becke */
#define XC_HYB_GGA_XC_B3LYP   402 /* The (in)famous B3LYP                  */
#define XC_HYB_GGA_XC_B3P86   403 /* Perdew 86 hybrid similar to B3PW91    */
#define XC_HYB_GGA_XC_mPW3PW  415 /* mixture with the mPW functional */
#define XC_HYB_GGA_XC_mPW3LYP 419 /* mixture of mPW and LYP */

/*************************************************************/
static void
gga_xc_b3_init(void *p_)
{
  const struct { FLOAT a0, ax, ac; int fldac, fggax, fggac; } 
  par[] = {
    {0.20,  0.72,  0.81,  XC_LDA_C_PW,      XC_GGA_X_B88,   XC_GGA_C_PW91}, /* b3pw91  */
    {0.20,  0.72,  0.81,  XC_LDA_C_VWN_RPA, XC_GGA_X_B88,   XC_GGA_C_LYP},  /* b3lyp   */
    {0.20,  0.72,  0.81,  XC_LDA_C_VWN_RPA, XC_GGA_X_B88,   XC_GGA_C_P86},  /* b3p86   */
    {0.20,  0.72,  0.81,  XC_LDA_C_VWN_RPA, XC_GGA_X_mPW91, XC_GGA_C_PW91}, /* mPW3PW  */
    {0.218, 0.709, 0.871, XC_LDA_C_VWN_RPA, XC_GGA_X_mPW91, XC_GGA_C_LYP}   /* mPW3LYP */
  };
  int func;

  XC(hyb_gga_type) *p = (XC(hyb_gga_type) *)p_;

  switch(p->info->number){
  case XC_HYB_GGA_XC_B3PW91:  func = 0; break;
  case XC_HYB_GGA_XC_B3LYP:   func = 1; break;
  case XC_HYB_GGA_XC_B3P86:   func = 2; break;
  case XC_HYB_GGA_XC_mPW3PW:  func = 3; break;
  case XC_HYB_GGA_XC_mPW3LYP: func = 4; break;
  default:
    fprintf(stderr, "Internal error in gga_xc_b3_init\n");
    exit(1);
    break;
  }

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 2;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  p->exx_coef = par[func].a0;

  XC(lda_init)(&p->mix->lda_mix[0], XC_LDA_X, p->nspin);
  p->mix->lda_coef[0] = 1.0 - par[func].a0 - par[func].ax;

  XC(lda_init)  (&p->mix->lda_mix[1], par[func].fldac, p->nspin);
  p->mix->lda_coef[1] = 1.0 - par[func].ac;

  /* set the vwn part with the spin interpolation scheme originally used in Gaussian */
  if(par[func].fldac == XC_LDA_C_VWN_RPA) /* b3lyp like functionals */
    XC(lda_c_vwn_set_params)(&p->mix->lda_mix[1], 1);

  XC(gga_init)(&p->mix->gga_mix[0], par[func].fggax, p->nspin);
  p->mix->gga_coef[0] = par[func].ax;
  XC(gga_init)(&p->mix->gga_mix[1], par[func].fggac, p->nspin);
  p->mix->gga_coef[1] = par[func].ac;
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3pw91) = {
  XC_HYB_GGA_XC_B3PW91,
  XC_EXCHANGE_CORRELATION,
  "B3PW91",
  XC_FAMILY_HYB_GGA,
  "AD Becke, J. Chem. Phys. 98, 5648 (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3lyp) = {
  XC_HYB_GGA_XC_B3LYP,
  XC_EXCHANGE_CORRELATION,
  "B3LYP",
  XC_FAMILY_HYB_GGA,
  "PJ Stephens, FJ Devlin, CF Chabalowski, MJ Frisch, J. Phys. Chem. 98 11623 (1994)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b3p86) = {
  XC_HYB_GGA_XC_B3P86,
  XC_EXCHANGE_CORRELATION,
  "B3P86",
  XC_FAMILY_HYB_GGA,
  "Defined through Gaussian implementation",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3pw) = {
  XC_HYB_GGA_XC_mPW3PW,
  XC_EXCHANGE_CORRELATION,
  "mPW3PW of Adamo & Barone",
  XC_FAMILY_GGA,
  "C Adamo and V Barone, J. Chem. Phys. 108, 664 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  gga_xc_b3_init, 
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_mpw3lyp) = {
  XC_HYB_GGA_XC_mPW3LYP,
  XC_EXCHANGE_CORRELATION,
  "mPW3LYP",
  XC_FAMILY_GGA,
  "Y Zhao and DGJ Truhlar, Phys. Chem. A 108, 6908 (2004)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_b3_init, 
  NULL, NULL, NULL
};
