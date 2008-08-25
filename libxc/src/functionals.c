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

#include <stdlib.h>

#include "xc.h"

extern XC(func_info_type) 
  *XC(lda_known_funct)[], 
  *XC(gga_known_funct)[],
  *XC(hyb_gga_known_funct)[],
  *XC(mgga_known_funct)[];

int XC(family_from_id)(int id)
{
  int i;

  /* first let us check if it is an LDA */
  for(i=0; XC(lda_known_funct)[i]!=NULL; i++){
    if(XC(lda_known_funct)[i]->number == id) return XC_FAMILY_LDA;
  }

  /* or is it a GGA? */
  for(i=0; XC(gga_known_funct)[i]!=NULL; i++){
    if(XC(gga_known_funct)[i]->number == id) return XC_FAMILY_GGA;
  }

  /* or is it a hybrid GGA? */
  for(i=0; XC(hyb_gga_known_funct)[i]!=NULL; i++){
    if(XC(hyb_gga_known_funct)[i]->number == id) return XC_FAMILY_HYB_GGA;
  }

  /* or is it a meta GGA? */
  for(i=0; XC(mgga_known_funct)[i]!=NULL; i++){
    if(XC(mgga_known_funct)[i]->number == id) return XC_FAMILY_MGGA;
  }

  return XC_FAMILY_UNKNOWN;
}
