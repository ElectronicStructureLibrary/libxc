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


#include "util.h"

#define XC_MGGA_C_CS          72 /* Colle and Salvetti */

/*
    [1] Eq. (15) in http://dx.doi.org/10.1103/PhysRevB.37.785
    [2] CS2 in http://www.molpro.net/info/2012.1/doc/manual/node192.html

  there is a gamma(r) in [1] absent in [2]. This should be irrelevant
  for spin unpolarized. In any case, it seems that even in that case,
  libxc does not give the same as molpro, but I am unable to
  understand why...
*/

#include "maple2c/mgga_c_cs.c"

#define func maple2c_func
#include "work_mgga_c.c"

const XC(func_info_type) XC(func_info_mgga_c_cs) = {
  XC_MGGA_C_CS,
  XC_CORRELATION,
  "Colle and Salvetti",
  XC_FAMILY_MGGA,
  {&xc_ref_Colle1975_329, &xc_ref_Lee1988_785, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  NULL, NULL, 
  NULL, NULL, work_mgga_c,
};
