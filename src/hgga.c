/*
 Copyright (C) 2006-2021 M.A.L. Marques
               2021 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"
#include "funcs_hgga.c"

void xc_evaluate_hgga(const xc_func_type *func, int max_order,
                const xc_input_variables *in, xc_output_variables *out)
{
  int ii, check;
  int orders[XC_MAXIMUM_ORDER+1] =
    {out->zk != NULL, out->vrho != NULL, out->v2rho2 != NULL,
     out->v3rho3 != NULL, out->v4rho4 != NULL};

  /* turn off orders smaller than max_order */
  for(ii=max_order+1; ii <= XC_MAXIMUM_ORDER; ii++)
    orders[ii] = 0;

  /* check if all variables make sense */
  check = xc_input_variables_sanity_check(in, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    fprintf(stderr, "Field %s is not allocated\n", xc_input_variables_name[check]);
    exit(1);
  }
  
  check = xc_output_variables_sanity_check(out, orders, func->info->family, func->info->flags);
  if(check >= 0){ /* error */
    if(check >= 1000)
      fprintf(stderr, "Functional does not provide an implementation of the %d-th derivative\n", check-1000);
    else
      fprintf(stderr, "Field %s is not allocated\n", xc_output_variables_name[check]);
    exit(1);
  }
  
  xc_output_variables_initialize(out, in->np, func->nspin);

  /* call the hGGA routines */
  if(func->info->hgga != NULL){
    if(func->nspin == XC_UNPOLARIZED){
      if(func->info->hgga->unpol[max_order] != NULL)
        func->info->hgga->unpol[max_order](func, in->np, in->rho, in->sigma, in->lapl, in->tau, in->exx, out);
    }else{
      if(func->info->hgga->pol[max_order] != NULL)
        func->info->hgga->pol[max_order](func, in->np, in->rho, in->sigma, in->lapl, in->tau, in->exx, out);
    }
  }

  if(func->mix_coef != NULL)
    xc_mix_func(func, in, out);
}
