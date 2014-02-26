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


/* this function converts the spin-density into total density and
	 relative magnetization */
/* inline */ void
XC(rho2dzeta)(int nspin, const FLOAT *rho, FLOAT *d, FLOAT *zeta)
{
  if(nspin==XC_UNPOLARIZED){
    *d    = max(rho[0], 0.0);
    *zeta = 0.0;
  }else{
    *d = rho[0] + rho[1];
    if(*d > 0.0){
      *zeta = (rho[0] - rho[1])/(*d);
      *zeta = min(*zeta,  1.0);
      *zeta = max(*zeta, -1.0);
    }else{
      *d    = 0.0;
      *zeta = 0.0;
    }
  }
}

/* inline */ void
XC(fast_fzeta)(const FLOAT x, const int nspin, const int order, FLOAT * fz){

  FLOAT aa, bb, aa2, bb2;

  if(nspin != XC_UNPOLARIZED){
    aa = CBRT(1.0 + x);
    bb = CBRT(1.0 - x);
    
    aa2 = aa*aa;
    bb2 = bb*bb;
    
    fz[0] = (aa2*aa2 + bb2*bb2 - 2.0)/FZETAFACTOR;
    if(order < 1) return;
    fz[1] = (aa - bb)*(4.0/3.0)/FZETAFACTOR;
    if(order < 2) return;
    fz[2] = ((4.0/9.0)/FZETAFACTOR)*(ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -2.0/3.0) + pow(1.0 - (x), -2.0/3.0)));
    if(order < 3) return;
    fz[3] = (-(8.0/27.0)/FZETAFACTOR)*(ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -5.0/3.0) - pow(1.0 - (x), -5.0/3.0)));
  } else {
    fz[0] = 0.0;
    fz[1] = 0.0;
    fz[2] = (8.0/9.0)/FZETAFACTOR;
    fz[3] = 0.0;
  }
}

/* initializes the mixing */
void 
XC(mix_init)(XC(func_type) *p, int n_funcs, const int *funcs_id, const FLOAT *mix_coef)
{
  int ii;

  assert(p != NULL);
  assert(p->func_aux == NULL && p->mix_coef == NULL);

  /* allocate structures needed for */
  p->n_func_aux = n_funcs;
  p->mix_coef   = (FLOAT *) malloc(n_funcs*sizeof(FLOAT));
  p->func_aux   = (XC(func_type) **) malloc(n_funcs*sizeof(XC(func_type) *));

  for(ii=0; ii<n_funcs; ii++){
    p->mix_coef[ii] = mix_coef[ii];
    p->func_aux[ii] = (XC(func_type) *) malloc(sizeof(XC(func_type)));
    XC(func_init) (p->func_aux[ii], funcs_id[ii], p->nspin);
  }
  
}

xc_gga_enhancement_t
XC(get_gga_enhancement_factor)(int func_id)
{
  switch(func_id){

  case XC_GGA_X_WC:
    return XC(gga_x_wc_enhance);

  case XC_GGA_X_PBE:
  case XC_GGA_X_PBE_R:
  case XC_GGA_X_PBE_SOL:
  case XC_GGA_X_XPBE:
  case XC_GGA_X_PBE_JSJR:
  case XC_GGA_X_PBEK1_VDW:
  case XC_GGA_X_RGE2:
  case XC_GGA_X_APBE:
  case XC_GGA_X_PBEINT:
  case XC_GGA_X_PBE_TCA:
    return XC(gga_x_pbe_enhance);

  case XC_GGA_X_PW91:
  case XC_GGA_X_MPW91:
    return XC(gga_x_pw91_enhance);

  case XC_GGA_X_RPBE:
    return XC(gga_x_rpbe_enhance);

  case XC_GGA_X_HTBS:
    return XC(gga_x_htbs_enhance);

  case XC_GGA_X_B86:
    return XC(gga_x_b86_enhance);

  case XC_GGA_X_B86_MGC:
    return XC(gga_x_b86_mgc_enhance);

  case XC_GGA_X_B88:
  case XC_GGA_X_OPTB88_VDW:
  case XC_GGA_X_MB88:
    return XC(gga_x_b88_enhance);

  case XC_GGA_X_G96:
    return XC(gga_x_g96_enhance);

  case XC_GGA_X_PW86:
  case XC_GGA_X_RPW86:
    return XC(gga_x_pw86_enhance);

  case XC_GGA_X_AIRY:
  case XC_GGA_X_LAG:
    return XC(gga_x_airy_enhance);

  case XC_GGA_X_BAYESIAN:
    return XC(gga_x_bayesian_enhance);

  case XC_GGA_X_BPCCAC:
    return XC(gga_x_bpccac_enhance);

  case XC_GGA_X_C09X:
    return XC(gga_x_c09x_enhance);

  case XC_GGA_X_AM05:
    return XC(gga_x_am05_enhance);

  case XC_GGA_X_DK87_R1:
  case XC_GGA_X_DK87_R2:
    return XC(gga_x_dk87_enhance);

  case XC_GGA_X_HERMAN:
    return XC(gga_x_herman_enhance);

  case XC_GGA_X_LG93:
    return XC(gga_x_lg93_enhance);

  case XC_GGA_X_LV_RPW86:
    return XC(gga_x_lv_rpw86_enhance);

  case XC_GGA_X_MPBE:
    return XC(gga_x_mpbe_enhance);

  case XC_GGA_X_OPTX:
    return XC(gga_x_optx_enhance);

  case XC_GGA_X_SOGGA11:
  case XC_HYB_GGA_X_SOGGA11_X:
    return XC(gga_x_sogga11_enhance);

  case XC_GGA_X_SSB_SW:
  case XC_GGA_X_SSB:
  case XC_GGA_X_SSB_D:
    return XC(gga_x_ssb_sw_enhance);

  case XC_GGA_X_VMT_PBE:
  case XC_GGA_X_VMT_GE:
  case XC_GGA_X_VMT84_PBE:
  case XC_GGA_X_VMT84_GE:
    return XC(gga_x_vmt_enhance);

  default:
    fprintf(stderr, "Internal error in get_gga_enhancement\n");
    exit(1);
  }
}
