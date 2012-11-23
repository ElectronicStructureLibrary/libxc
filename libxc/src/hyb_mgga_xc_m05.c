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

#define XC_HYB_MGGA_XC_M05      438 /* M05 functional of Minnesota                      */
#define XC_HYB_MGGA_XC_M05_2X   439 /* M05-2X functional of Minnesota                   */
#define XC_HYB_MGGA_XC_BX88B95  440 /* Mixture of B88 with BC95 (B1B95)                 */
#define XC_HYB_MGGA_XC_BX86B95  441 /* Mixture of B86 with BC95                         */
#define XC_HYB_MGGA_XC_PWX86B95 442 /* Mixture of PW86 with BC95                        */
#define XC_HYB_MGGA_XC_BB1K     443 /* Mixture of B88 with BC95 from Zhao and Truhlar   */
#define XC_HYB_MGGA_XC_MPW1B95  445 /* Mixture of mPW91 with BC95 from Zhao and Truhlar */
#define XC_HYB_MGGA_XC_MPWB1K   446 /* Mixture of mPW91 with BC95 for kinetics          */
#define XC_HYB_MGGA_XC_X1B95    447 /* Mixture of X with BC95                           */
#define XC_HYB_MGGA_XC_XB1K     448 /* Mixture of X with BC95 for kinetics              */
#define XC_HYB_MGGA_XC_VSXC     444 /* */

/*************************************************************/
void
XC(hyb_mgga_xc_m05_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_M05, XC_MGGA_C_M05};
  static FLOAT funcs_coef[2] = {1.0 - 0.28, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.28;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_m05) = {
  XC_HYB_MGGA_XC_M05,
  XC_EXCHANGE_CORRELATION,
  "M05 functional of Minnesota",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao, NE Schultz, and DG Truhlar, J. Chem. Phys. 123, 161103 (2005)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_m05_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_m05_2x_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_MGGA_X_M05_2X, XC_MGGA_C_M05_2X};
  static FLOAT funcs_coef[2] = {1.0 - 0.56, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.56;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_m05_2x) = {
  XC_HYB_MGGA_XC_M05_2X,
  XC_EXCHANGE_CORRELATION,
  "M05-2X functional of Minnesota",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao, NE Schultz, and DG Truhlar, J. Chem. Theory Comput. 2, 364 (2006)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_m05_2x_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_bx88bc95_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.28, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.28;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_bx88bc95) = {
  XC_HYB_MGGA_XC_BX88B95,
  XC_EXCHANGE_CORRELATION,
  "Mixture of B88 with BC95 (B1B95)",
  XC_FAMILY_HYB_MGGA,
  "A Becke, J. Chem. Phys. 104, 1040 (1996)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_bx88bc95_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_bx86bc95_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B86, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.28, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.28;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_bx86bc95) = {
  XC_HYB_MGGA_XC_BX86B95,
  XC_EXCHANGE_CORRELATION,
  "Mixture of B86 with BC95",
  XC_FAMILY_HYB_MGGA,
  "A Becke, J. Chem. Phys. 104, 1040 (1996)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_bx86bc95_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_pwx86bc95_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_PW86, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.29, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.29;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_pwx86bc95) = {
  XC_HYB_MGGA_XC_PWX86B95,
  XC_EXCHANGE_CORRELATION,
  "Mixture of PW86 with BC95",
  XC_FAMILY_HYB_MGGA,
  "A Becke, J. Chem. Phys. 104, 1040 (1996)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_pwx86bc95_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_bb1k_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_B88, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.42, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.42;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_bb1k) = {
  XC_HYB_MGGA_XC_BB1K,
  XC_EXCHANGE_CORRELATION,
  "Mixture of B88 with BC95 from Zhao and Truhlar",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao and DG Truhlar, J. Phys. Chem. A 108, 2715 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_bb1k_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_mpw1b95_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_MPW91, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.31, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.31;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_mpw1b95) = {
  XC_HYB_MGGA_XC_MPW1B95,
  XC_EXCHANGE_CORRELATION,
  "Mixture of mPW91 with BC95 from Zhao and Truhlar",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao and DG Truhlar, J. Phys. Chem. A 108, 6908-6918 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_mpw1b95_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_mpwb1k_init)(XC(func_type) *p)
{
  static int   funcs_id  [2] = {XC_GGA_X_MPW91, XC_MGGA_C_BC95};
  static FLOAT funcs_coef[2] = {1.0 - 0.44, 1.0};

  XC(mix_init)(p, 2, funcs_id, funcs_coef);
  p->cam_alpha = 0.44;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_mpwb1k) = {
  XC_HYB_MGGA_XC_MPWB1K,
  XC_EXCHANGE_CORRELATION,
  "Mixture of mPW91 with BC95 for kinetics",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao and DG Truhlar, J. Phys. Chem. A 108, 6908-6918 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_mpwb1k_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_x1b95_init)(XC(func_type) *p)
{
  const FLOAT a1=0.675, a2=0.235, a0=0.30;

  static int   funcs_id  [3] = {XC_GGA_X_B88, XC_GGA_X_PW91, XC_MGGA_C_BC95};
  FLOAT funcs_coef[3];

  funcs_coef[0] = a0*a1;
  funcs_coef[1] = a0*a2;
  funcs_coef[2] = 1.0;

  XC(mix_init)(p, 3, funcs_id, funcs_coef);
  p->cam_alpha = a0;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_x1b95) = {
  XC_HYB_MGGA_XC_X1B95,
  XC_EXCHANGE_CORRELATION,
  "Mixture of X with BC95",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao and DG Truhlar, J. Phys. Chem. A 108, 6908-6918 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_x1b95_init),
  NULL, NULL, NULL, NULL,
};


/*************************************************************/
void
XC(hyb_mgga_xc_xb1k_init)(XC(func_type) *p)
{
  const FLOAT a1=0.675, a2=0.235, a0=0.43;

  static int   funcs_id  [3] = {XC_GGA_X_B88, XC_GGA_X_PW91, XC_MGGA_C_BC95};
  FLOAT funcs_coef[3];

  funcs_coef[0] = a0*a1;
  funcs_coef[1] = a0*a2;
  funcs_coef[2] = 1.0;

  XC(mix_init)(p, 3, funcs_id, funcs_coef);
  p->cam_alpha = a0;
}

XC(func_info_type) XC(func_info_hyb_mgga_xc_xb1k) = {
  XC_HYB_MGGA_XC_XB1K,
  XC_EXCHANGE_CORRELATION,
  "Mixture of X with BC95 for kinetics",
  XC_FAMILY_HYB_MGGA,
  "Y Zhao and DG Truhlar, J. Phys. Chem. A 108, 6908-6918 (2004)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  XC(hyb_mgga_xc_xb1k_init),
  NULL, NULL, NULL, NULL,
};
