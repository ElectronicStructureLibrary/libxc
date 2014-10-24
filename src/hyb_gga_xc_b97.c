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

#define XC_HYB_GGA_XC_B97      407 /* Becke 97                                 */
#define XC_HYB_GGA_XC_B97_1    408 /* Becke 97-1                               */
#define XC_HYB_GGA_XC_B97_2    410 /* Becke 97-2                               */
#define XC_HYB_GGA_XC_B97_K    413 /* Boese-Martin for Kinetics                */
#define XC_HYB_GGA_XC_B97_3    414 /* Becke 97-3                               */
#define XC_HYB_GGA_XC_SB98_1a  420 /* Schmider-Becke 98 parameterization 1a    */
#define XC_HYB_GGA_XC_SB98_1b  421 /* Schmider-Becke 98 parameterization 1b    */
#define XC_HYB_GGA_XC_SB98_1c  422 /* Schmider-Becke 98 parameterization 1c    */
#define XC_HYB_GGA_XC_SB98_2a  423 /* Schmider-Becke 98 parameterization 2a    */
#define XC_HYB_GGA_XC_SB98_2b  424 /* Schmider-Becke 98 parameterization 2b    */
#define XC_HYB_GGA_XC_SB98_2c  425 /* Schmider-Becke 98 parameterization 2c    */
#define XC_HYB_GGA_XC_WB97     463 /* Chai and Head-Gordon                     */
#define XC_HYB_GGA_XC_WB97X    464 /* Chai and Head-Gordon                     */

static void
hyb_gga_xc_b97_init(XC(func_type) *p)
{
  const struct { int iGGA; FLOAT a0; } par[] = {
    {XC_GGA_XC_B97,     0.1943},
    {XC_GGA_XC_B97_1,   0.21},
    {XC_GGA_XC_B97_2,   0.21},
    {XC_GGA_XC_B97_K,   0.42},
    {XC_GGA_XC_B97_3,   2.692880E-01},
    {XC_GGA_XC_SB98_1a, 0.229015},
    {XC_GGA_XC_SB98_1b, 0.199352},
    {XC_GGA_XC_SB98_1c, 0.192416},
    {XC_GGA_XC_SB98_2a, 0.232055},
    {XC_GGA_XC_SB98_2b, 0.237978},
    {XC_GGA_XC_SB98_2c, 0.219847},
    {XC_GGA_XC_WB97,    1.000000},
    {XC_GGA_XC_WB97X,   1.000000}
  };
  
  int func;
  FLOAT one = 1.0;
  
  switch(p->info->number){
  case XC_HYB_GGA_XC_B97:      func =  0; break;
  case XC_HYB_GGA_XC_B97_1:    func =  1; break;
  case XC_HYB_GGA_XC_B97_2:    func =  2; break;
  case XC_HYB_GGA_XC_B97_K:    func =  3; break;
  case XC_HYB_GGA_XC_B97_3:    func =  4; break;
  case XC_HYB_GGA_XC_SB98_1a:  func =  5; break;
  case XC_HYB_GGA_XC_SB98_1b:  func =  6; break;
  case XC_HYB_GGA_XC_SB98_1c:  func =  7; break;
  case XC_HYB_GGA_XC_SB98_2a:  func =  8; break;
  case XC_HYB_GGA_XC_SB98_2b:  func =  9; break;
  case XC_HYB_GGA_XC_SB98_2c:  func = 10; break;
  case XC_HYB_GGA_XC_WB97:     func = 11; break;
  case XC_HYB_GGA_XC_WB97X:    func = 12; break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_xc_b97_init\n");
    exit(1);
    break;
  }

  XC(mix_init)(p, 1, &(par[func].iGGA), &one);
  p->cam_alpha = par[func].a0;

  if(p->info->number == XC_HYB_GGA_XC_WB97) {
    p->cam_omega =  0.4;
    p->cam_beta  = -1.0;
  }else if(p->info->number == XC_HYB_GGA_XC_WB97X) {
    p->cam_omega =  0.3;
    p->cam_beta  = -(1.0 - 1.57706e-01);
  }
}

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97) = {
  XC_HYB_GGA_XC_B97,
  XC_EXCHANGE_CORRELATION,
  "Becke 97",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Becke1997_8554, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init,
  NULL, NULL, NULL, NULL /* this is taken care by the generic routine */
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_1) = {
  XC_HYB_GGA_XC_B97_1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-1",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init,
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_2) = {
  XC_HYB_GGA_XC_B97_2,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Wilson2001_9233, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_k) = {
  XC_HYB_GGA_XC_B97_K,
  XC_EXCHANGE_CORRELATION,
  "Boese-Martin for Kinetics",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_b97_3) = {
  XC_HYB_GGA_XC_B97_3,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-3",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Keal2005_121103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1a) = {
  XC_HYB_GGA_XC_SB98_1a,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1b) = {
  XC_HYB_GGA_XC_SB98_1b,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_1c) = {
  XC_HYB_GGA_XC_SB98_1c,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2a) = {
  XC_HYB_GGA_XC_SB98_2a,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2b) = {
  XC_HYB_GGA_XC_SB98_2b,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_sb98_2c) = {
  XC_HYB_GGA_XC_SB98_2c,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_wb97) = {
  XC_HYB_GGA_XC_WB97,
  XC_EXCHANGE_CORRELATION,
  "wB97 range-separated functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Chai2008_084106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};

const XC(func_info_type) XC(func_info_hyb_gga_xc_wb97x) = {
  XC_HYB_GGA_XC_WB97X,
  XC_EXCHANGE_CORRELATION,
  "wB97X range-separated functional",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Chai2008_084106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_xc_b97_init, 
  NULL, NULL, NULL, NULL
};
