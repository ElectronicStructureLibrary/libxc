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

extern xc_func_info_type 
  *lda_known_funct[], 
  *gga_known_funct[],
  *hyb_gga_known_funct[];

int xc_family_from_id(int id)
{
  int i;

  /* first let us check if it is an LDA */
  for(i=0; lda_known_funct[i]!=NULL; i++){
    if(lda_known_funct[i]->number == id) return XC_FAMILY_LDA;
  }

  /* or is it a GGA? */
  for(i=0; gga_known_funct[i]!=NULL; i++){
    if(gga_known_funct[i]->number == id) return XC_FAMILY_GGA;
  }

  /* or is it a hybrid GGA? */
  for(i=0; hyb_gga_known_funct[i]!=NULL; i++){
    if(hyb_gga_known_funct[i]->number == id) return XC_FAMILY_HYB_GGA;
  }

  return XC_FAMILY_UNKNOWN;
}
