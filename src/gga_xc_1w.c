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

#define XC_GGA_XC_PBE1W    173  /* Functionals fitted for water */
#define XC_GGA_XC_MPWLYP1W 174  /* Functionals fitted for water */
#define XC_GGA_XC_PBELYP1W 175  /* Functionals fitted for water */

static void
gga_xc_pbe1w_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 1;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  XC(lda_init)(&p->mix->lda_mix[0], XC_LDA_C_VWN, p->nspin);
  p->mix->lda_coef[0] = (1.0 - 74.0/100.0);

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_PBE, p->nspin);
  p->mix->gga_coef[0] = 1.0;

  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_C_PBE, p->nspin);
  p->mix->gga_coef[1] = 74.0/100.0;
}

const XC(func_info_type) XC(func_info_gga_xc_pbe1w) = {
  XC_GGA_XC_PBE1W,
  XC_EXCHANGE_CORRELATION,
  "PBE1W",
  XC_FAMILY_GGA,
  "EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  gga_xc_pbe1w_init, 
  NULL, NULL, NULL
};

static void
gga_xc_mpwlyp1w_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 1;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  XC(lda_init)(&p->mix->lda_mix[0], XC_LDA_C_VWN, p->nspin);
  p->mix->lda_coef[0] = (1.0 - 88.0/100.0);

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_mPW91, p->nspin);
  p->mix->gga_coef[0] = 1.0;

  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_C_LYP, p->nspin);
  p->mix->gga_coef[1] = 88.0/100.0;
}

const XC(func_info_type) XC(func_info_gga_xc_mpwlyp1w) = {
  XC_GGA_XC_MPWLYP1W,
  XC_EXCHANGE_CORRELATION,
  "mPWLYP1w",
  XC_FAMILY_GGA,
  "EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_mpwlyp1w_init, 
  NULL, NULL, NULL
};

static void
gga_xc_pbelyp1w_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->mix = (XC(mix_func_type) *) malloc(sizeof(XC(mix_func_type)));
  XC(mix_func_init)(p->mix, XC_FAMILY_GGA, p->nspin);

  p->mix->lda_n = 1;
  p->mix->gga_n = 2;
  XC(mix_func_alloc)(p->mix);

  XC(lda_init)(&p->mix->lda_mix[0], XC_LDA_C_VWN, p->nspin);
  p->mix->lda_coef[0] = (1.0 - 74.0/100.0);

  XC(gga_init)(&p->mix->gga_mix[0], XC_GGA_X_PBE, p->nspin);
  p->mix->gga_coef[0] = 1.0;

  XC(gga_init)(&p->mix->gga_mix[1], XC_GGA_C_LYP, p->nspin);
  p->mix->gga_coef[1] = 74.0/100.0;
}

const XC(func_info_type) XC(func_info_gga_xc_pbelyp1w) = {
  XC_GGA_XC_PBELYP1W,
  XC_EXCHANGE_CORRELATION,
  "PBELYP1W",
  XC_FAMILY_GGA,
  "EE Dahlke and DG Truhlar, J. Phys. Chem. B 109, 15677 (2005)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_pbelyp1w_init, 
  NULL, NULL, NULL
};

