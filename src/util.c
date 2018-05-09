/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"


/* this function converts the spin-density into total density and
	 relative magnetization */
/* inline */ void
xc_rho2dzeta(int nspin, const double *rho, double *d, double *zeta)
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

xc_gga_enhancement_t
xc_get_gga_enhancement_factor(int func_id)
{
  switch(func_id){

  case XC_GGA_X_WC:
    return xc_gga_x_wc_enhance;

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
    return xc_gga_x_pbe_enhance;

  case XC_GGA_X_PW91:
  case XC_GGA_X_MPW91:
    return xc_gga_x_pw91_enhance;

  case XC_GGA_X_RPBE:
    return xc_gga_x_rpbe_enhance;

  case XC_GGA_X_HTBS:
    return xc_gga_x_htbs_enhance;

  case XC_GGA_X_B86:
  case XC_GGA_X_B86_MGC:
  case XC_GGA_X_B86_R:
    return xc_gga_x_b86_enhance;

  case XC_GGA_X_B88:
  case XC_GGA_X_OPTB88_VDW:
  case XC_GGA_X_MB88:
    return xc_gga_x_b88_enhance;

  case XC_GGA_X_G96:
    return xc_gga_x_g96_enhance;

  case XC_GGA_X_PW86:
  case XC_GGA_X_RPW86:
    return xc_gga_x_pw86_enhance;

  case XC_GGA_X_AIRY:
  case XC_GGA_X_LAG:
    return xc_gga_x_airy_enhance;

  case XC_GGA_X_BAYESIAN:
    return xc_gga_x_bayesian_enhance;

  case XC_GGA_X_BPCCAC:
    return xc_gga_x_bpccac_enhance;

  case XC_GGA_X_C09X:
    return xc_gga_x_c09x_enhance;

  case XC_GGA_X_AM05:
    return xc_gga_x_am05_enhance;

  case XC_GGA_X_DK87_R1:
  case XC_GGA_X_DK87_R2:
    return xc_gga_x_dk87_enhance;

  case XC_GGA_X_HERMAN:
    return xc_gga_x_herman_enhance;

  case XC_GGA_X_LG93:
    return xc_gga_x_lg93_enhance;

  case XC_GGA_X_LV_RPW86:
    return xc_gga_x_lv_rpw86_enhance;

  case XC_GGA_X_MPBE:
    return xc_gga_x_mpbe_enhance;

  case XC_GGA_X_OPTX:
    return xc_gga_x_optx_enhance;

  case XC_GGA_X_SOGGA11:
  case XC_HYB_GGA_X_SOGGA11_X:
    return xc_gga_x_sogga11_enhance;

  case XC_GGA_X_SSB_SW:
  case XC_GGA_X_SSB:
  case XC_GGA_X_SSB_D:
    return xc_gga_x_ssb_sw_enhance;

  case XC_GGA_X_VMT_PBE:
  case XC_GGA_X_VMT_GE:
  case XC_GGA_X_VMT84_PBE:
  case XC_GGA_X_VMT84_GE:
    return xc_gga_x_vmt_enhance;

  default:
    fprintf(stderr, "Internal error in get_gga_enhancement\n");
    exit(1);
  }
}


const char *get_kind(const xc_func_type *func) {
  switch(func->info->kind) {
    case(XC_EXCHANGE):
      return "XC_EXCHANGE";

    case(XC_CORRELATION):
      return "XC_CORRELATION";

    case(XC_EXCHANGE_CORRELATION):
      return "XC_EXCHANGE_CORRELATION";

    case(XC_KINETIC):
      return "XC_KINETIC";

    default:
      printf("Internal error in get_kind.\n");
      return "";
  }
}

const char *get_family(const xc_func_type *func) {
  switch(func->info->family) {
    case(XC_FAMILY_UNKNOWN):
      return "XC_FAMILY_UNKNOWN";

    case(XC_FAMILY_LDA):
      return "XC_FAMILY_LDA";

    case(XC_FAMILY_GGA):
      return "XC_FAMILY_GGA";

    case(XC_FAMILY_MGGA):
      return "XC_FAMILY_MGGA";

    case(XC_FAMILY_LCA):
      return "XC_FAMILY_LCA";

    case(XC_FAMILY_OEP):
      return "XC_FAMILY_OEP";

    case(XC_FAMILY_HYB_GGA):
      return "XC_FAMILY_HYB_GGA";

    case(XC_FAMILY_HYB_MGGA):
      return "XC_FAMILY_HYB_MGGA";

    default:
      printf("Internal error in get_family.\n");
      return "";
  }
}

/* this function checks if it should use the default or
   the user assigned value for an external parameter */
double
get_ext_param(const func_params_type *params, const double *values, int index)
{
  FILE *par_in;
  int ii, nn;
  double dd;
  
  /* 
     If libxc finds a file in the current directory name
     "libxc.params", it will try to read the parameters for the
     current functional from it. This file should contain one
     parameter per line. E.g., for the x_pbe functional:

       ------------------ <start libxc.params>
       0.8040              # _kappa
       0.2195149727645171  # _mu (PBE)
       ------------------ <end libxc.params>

     Note that this only works for functionals whose parameters can be
     set by set_ext_params.
  */
  /* Commented as considered dangerous ;)
  if((par_in = fopen("libxc.params","rb"))){
    for(ii=0; ii<index; ii++)
      fscanf(par_in, "%*[^\n]\n", NULL);

    nn = fscanf(par_in, "%lf", &dd);
    fclose(par_in);

    if(nn == 1)
      return dd;
  }
  */

  if(values == NULL || values[index] == XC_EXT_PARAMS_DEFAULT)
    return params[index].value; /* return default value */
  else
    return values[index]; /* return user assigned value */
}
