/*
 Copyright (C) 2017 Susi Lehtola

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

#define XC_HYB_MGGA_XC_MPW1KCIS    566 /* Modified Perdew-Wang + KCIS hybrid */
#define XC_HYB_MGGA_XC_MPWKCIS1K   567 /* Modified Perdew-Wang + KCIS hybrid with more exact exchange */
#define XC_HYB_MGGA_XC_PBE1KCIS    568 /* Perdew-Burke-Ernzerhof + KCIS hybrid */
#define XC_HYB_MGGA_XC_TPSS1KCIS   569 /* TPSS hybrid with KCIS correlation */

static void
hyb_mgga_xc_kcis_init(XC(func_type) *p)
{
  /* Exchange functional */
  int xid;
  /* Fraction of exact exchange */
  FLOAT exx;
  /* Array */
  int funcs_id[2];
  FLOAT funcs_coef[2];

  switch(p->info->number){
  case XC_HYB_MGGA_XC_MPW1KCIS:
    xid=XC_GGA_X_MPW91;
    exx=0.15;
    break;
  case XC_HYB_MGGA_XC_MPWKCIS1K:
    xid=XC_GGA_X_MPW91;
    exx=0.41;
    break;
  case XC_HYB_MGGA_XC_PBE1KCIS:
    xid=XC_GGA_X_PBE;
    exx=0.22;
    break;
  case XC_HYB_MGGA_XC_TPSS1KCIS:
    xid=XC_MGGA_X_TPSS;
    exx=0.13;
    break;
  default:
    fprintf(stderr, "Internal error in hyb_mgga_xc_kcis\n");
    exit(1);
  }

  /* Initialize mix */
  funcs_id[0] = xid;
  funcs_coef[0] = 1.0-exx;

  funcs_id[1] = XC_MGGA_C_KCIS;
  funcs_coef[1] = 1.0;
  
  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = exx;
}

const XC(func_info_type) XC(func_info_hyb_mgga_xc_mpw1kcis) = {
  XC_HYB_MGGA_XC_MPW1KCIS,
  XC_EXCHANGE_CORRELATION,
  "MPW1KCIS for barrier heights",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_2012, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_kcis_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

const XC(func_info_type) XC(func_info_hyb_mgga_xc_mpwkcis1k) = {
  XC_HYB_MGGA_XC_MPWKCIS1K,
  XC_EXCHANGE_CORRELATION,
  "MPWKCIS1K for barrier heights",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_2012, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_kcis_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

const XC(func_info_type) XC(func_info_hyb_mgga_xc_pbe1kcis) = {
  XC_HYB_MGGA_XC_PBE1KCIS,
  XC_EXCHANGE_CORRELATION,
  "PBE1KCIS for binding energies",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_415, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_kcis_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

const XC(func_info_type) XC(func_info_hyb_mgga_xc_tpss1kcis) = {
  XC_HYB_MGGA_XC_TPSS1KCIS,
  XC_EXCHANGE_CORRELATION,
  "TPSS1KCIS for thermochemistry and kinetics",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Zhao2005_43, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32,
  0, NULL, NULL,
  hyb_mgga_xc_kcis_init,
  NULL, NULL, NULL, NULL /* this is taken care of by the generic routine */
};

