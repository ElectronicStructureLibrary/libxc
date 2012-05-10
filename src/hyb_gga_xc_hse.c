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

#define XC_HYB_GGA_XC_HSE03 427 /* the 2003 version of the screened hybrid HSE */
#define XC_HYB_GGA_XC_HSE06 428 /* the 2006 version of the screened hybrid HSE */

static void
hyb_gga_xc_hse_init(void *p_)
{
  static int   funcs_id  [3] = {XC_GGA_X_WPBEH, XC_GGA_X_WPBEH, XC_GGA_C_PBE};
  static FLOAT funcs_coef[3] = {1.0, -0.25, 1.0};  
  XC(gga_type) *p = (XC(gga_type) *)p_;

  XC(gga_init_mix)(p, 3, funcs_id, funcs_coef);
  p->exx_coef = 0.25;
  
  /* Note that there is an enormous mess in the literature concerning
     the values of omega in HSE. This is due to an error in the
     original paper that stated that they had used omega=0.15. This
     was in fact not true, and the real value used was omega^HF =
     0.15/sqrt(2) ~ 0.1061 and omega^PBE = 0.15*cbrt(2) ~ 0.1890. In
     2006 Krukau et al [JCP 125, 224106 (2006)] tried to clarify the
     situation, and called HSE03 to the above choice of parameters,
     and called HSE06 to the functional where omega^HF=omega^PBE. By
     testing several properties for atoms they reached the conclusion
     that the best value for omega=0.11.

     Of course, codes are just as messy as the papers. In espresso
     HSE06 has the value omega=0.106. VASP, on the other hand, uses
     for HSE03 the same value omega^HF = omega^PBE = 0.3 (A^-1) ~
     0.1587 and for HSE06 omega^HF = omega^PBE = 0.2 (A^-1) ~ 0.1058.

     We try to follow the original definition of the functional. The
     default omega for XC_GGA_X_WPBEH is zero, so WPBEh reduces to
     PBEh
   */
  switch(p->info->number){
  case XC_HYB_GGA_XC_HSE03:
    /* in this case one should use omega^HF = 0.15/sqrt(2) and
       omega^PBE = 0.15*CBRT(2.0)*/
    XC(hyb_gga_xc_hse_set_params_)(p, 0.15*CBRT(2.0));
    break;
  case XC_HYB_GGA_XC_HSE06:
    /* in this case one should use omega^HF = omega^PBE = 0.11 */
    XC(hyb_gga_xc_hse_set_params_)(p, 0.11);
    break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_hse\n");
    exit(1);
  }

  p->exx_coef = 0.25;
}


void 
XC(hyb_gga_xc_hse_set_params)(XC(func_type) *p, FLOAT omega)
{
  assert(p != NULL && p->gga != NULL);
  XC(hyb_gga_xc_hse_set_params_)(p->gga, omega);
}


void 
XC(hyb_gga_xc_hse_set_params_)(XC(gga_type) *p, FLOAT omega)
{
  assert(p->func_aux[1] != NULL);
   (p->params);

   XC(gga_x_pbe_sr_set_params)(p->func_aux[1], omega);
}


const XC(func_info_type) XC(func_info_hyb_gga_xc_hse03) = {
  XC_HYB_GGA_XC_HSE03,
  XC_EXCHANGE_CORRELATION,
  "HSE03",
  XC_FAMILY_HYB_GGA,
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)\n"
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 124, 219906 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_hse06) = {
  XC_HYB_GGA_XC_HSE06,
  XC_EXCHANGE_CORRELATION,
  "HSE06",
  XC_FAMILY_HYB_GGA,
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)\n"
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 124, 219906 (2006)\n"
  "AV Krukau, OA Vydrov, AF Izmaylov, and GE Scuseria, J. Chem. Phys. 125, 224106 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_hse_init,
  NULL, NULL, NULL
};
