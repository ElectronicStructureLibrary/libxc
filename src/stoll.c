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

void 
XC(lda_stoll) (const XC(func_type) *pw, FLOAT dens, FLOAT zeta, int order, XC(lda_work_t) res[3])
{
  static const FLOAT sign[2] = {1.0, -1.0};
  int is;
  FLOAT opz[2] = {1.0 + zeta, 1.0 - zeta};

  res[2].rs[1]  = RS(dens);

  /* first we get the parallel contributions */
  for(is=0; is<2; is++){
    FLOAT opz13;

    if(opz[is] < pw->info->min_zeta){
      res[is].zk = 0.0;
      if(order >= 1) res[is].dedz = res[is].dedrs = 0.0;
    }else{
      FLOAT drssdrs, drssdz, d2rssdrsz, d2rssdz2;
      FLOAT LDA_zk, LDA_dedrs, LDA_d2edrs2;

      opz13 = CBRT(opz[is]);

      res[is].rs[1] = RS(dens*opz[is]/2.0);
      res[is].rs[0] = sqrt(res[is].rs[1]);
      res[is].rs[2] = res[is].rs[1]*res[is].rs[1];
      res[is].zeta  = sign[is];
      res[is].order = order;
  
      XC(lda_c_pw_func)(pw, &(res[is]));

      LDA_zk = res[is].zk;

      res[is].zk *= opz[is]/2.0;
      
      if(order < 1) continue;

      LDA_dedrs = res[is].dedrs;
      drssdrs   = M_CBRT2/opz13;
      drssdz    = -sign[is]*res[is].rs[1]/(3.0*opz[is]);

      res[is].dedrs = LDA_dedrs*drssdrs*opz[is]/2.0;
      res[is].dedz  = LDA_zk*sign[is]/2.0 + LDA_dedrs*drssdz*opz[is]/2.0;

      if(order < 2) continue;

      LDA_d2edrs2 = res[is].d2edrs2;
      d2rssdrsz   = -sign[is]*M_CBRT2/(3.0*opz13*opz[is]);
      d2rssdz2    = res[is].rs[1]*4.0/(9.0*opz[is]*opz[is]);

      res[is].d2edrs2 = LDA_d2edrs2*drssdrs*drssdrs*opz[is]/2.0;
      res[is].d2edrsz = sign[is]*LDA_dedrs*drssdrs/2.0 + (LDA_d2edrs2*drssdrs*drssdz + LDA_dedrs*d2rssdrsz)*opz[is]/2.0;
      res[is].d2edz2  = sign[is]*LDA_dedrs*drssdz + (LDA_d2edrs2*drssdz*drssdz + LDA_dedrs*d2rssdz2)*opz[is]/2.0;
    }
  }

  /* and now the perpendicular */
  res[2].rs[0]  = sqrt(res[2].rs[1]);
  res[2].rs[2]  = res[2].rs[1]*res[2].rs[1];
  res[is].zeta  = zeta;
  res[is].order = order;

  XC(lda_c_pw_func)(pw, &(res[2]));

  res[2].zk -= res[0].zk + res[1].zk;

  if(order < 1) return;

  res[2].dedrs -= res[0].dedrs + res[1].dedrs;
  res[2].dedz  -= res[0].dedz  + res[1].dedz;

  if(order < 2) return;

  res[2].d2edrs2 -= res[0].d2edrs2 + res[1].d2edrs2;
  res[2].d2edrsz -= res[0].d2edrsz + res[1].d2edrsz;
  res[2].d2edz2  -= res[0].d2edz2  + res[1].d2edz2;
}
