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

#define XC_MGGA_XC_CC06          229 /* Cancio and Chou 2006 */


static void 
mgga_xc_cc06_init(XC(func_type) *p)
{
  assert(p != NULL);

  p->n_func_aux  = 2;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));
  p->func_aux[1] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_X,    p->nspin);
  XC(func_init)(p->func_aux[1], XC_LDA_C_PW, p->nspin);
}


static void 
func(const XC(func_type) *pt, XC(mgga_work_c_t) *r)
{
  static FLOAT alpha = -0.0007, beta = 0.0080, gamma = 0.026;
  XC(lda_work_t) lda_x, lda_pw;
  FLOAT l_cnst, opz, omz, opz13, omz13, opz23, omz23, l, fxc_n, fxc_d, fxc;
  FLOAT dldz, dldus[2], dfxc;

  lda_pw.order = lda_x.order = r->order;
  lda_pw.rs    = lda_x.rs    = r->rs;
  lda_pw.zeta  = lda_x.zeta  = r->zeta;

  XC(lda_x_func)   (pt->func_aux[0], &lda_x);
  XC(lda_c_pw_func)(pt->func_aux[1], &lda_pw);

  l_cnst = CBRT(3.0/(2.0*4.0*M_PI));
  l_cnst = l_cnst*l_cnst/2.0;

  opz = 1.0 + r->zeta;
  omz = 1.0 - r->zeta;

  opz13 = CBRT(opz); opz23 = opz13*opz13;
  omz13 = CBRT(omz); omz23 = omz13*omz13;

  l     = l_cnst*(r->us[0]*opz*opz23 + r->us[1]*omz*omz23);

  fxc_n = alpha + beta*l;
  fxc_d = 1.0 + gamma*l;
  fxc   = 1.0 + fxc_n/fxc_d;

  r->f = (lda_x.e + lda_pw.e)*fxc;

  if(r->order < 1) return;

  dldz     = l_cnst*(5.0/3.0)*(r->us[0]*opz23 - r->us[1]*omz23);
  dldus[0] = l_cnst*opz*opz23;
  dldus[1] = l_cnst*omz*omz23;

  dfxc = -(alpha*gamma - beta)/(fxc_d*fxc_d);

  r->dfdrs    = (lda_x.dedrs + lda_pw.dedrs)*fxc;
  r->dfdz     = (lda_x.dedz + lda_pw.dedz)*fxc + (lda_x.e + lda_pw.e)*dfxc*dldz;
  r->dfdxt    = 0.0;
  r->dfdxs[0] = 0.0;
  r->dfdxs[1] = 0.0;
  r->dfdts[0] = 0.0;
  r->dfdts[1] = 0.0;
  r->dfdus[0] = (lda_x.e + lda_pw.e)*dfxc*dldus[0];
  r->dfdus[1] = (lda_x.e + lda_pw.e)*dfxc*dldus[1];

  if(r->order < 2) return;
}


#include "work_mgga_c.c"


const XC(func_info_type) XC(func_info_mgga_xc_cc06) = {
  XC_MGGA_XC_CC06,
  XC_EXCHANGE_CORRELATION,
  "Cancio and Chou 2006",
  XC_FAMILY_MGGA,
  {&xc_ref_Cancio2006_081202, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_LAPLACIAN | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  0, NULL, NULL,
  mgga_xc_cc06_init,
  NULL, NULL, NULL,
  work_mgga_c,
};

