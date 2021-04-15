/*
 Copyright (C) 2017 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xc.h"

int main(void) {
  int i, N, error;
  xc_func_type func;

  /* List of functionals */
  int *flist;
  char *p;

  /* Get list of available functionals */
  N = xc_number_of_functionals();
  flist = (int *) malloc(N*sizeof(int));
  xc_available_functional_numbers(flist);

  /* Loop over functionals */
  for(i=0;i<N;i++) {
    /* Functional id and name */
    int func_id;
    char *fname;
    /* Kind and family of functional */
    char kind[5], family[10];

    printf("Checking functional %i -> %i\n", i, flist[i]);

    func_id = flist[i];
    fname = xc_functional_get_name(func_id);

    /* Initialize functional */
    error = xc_func_init(&func, func_id, XC_UNPOLARIZED);
    if(error) {
      printf("Error initializing functional %i.\n", func_id);
      fflush(stdout);
      return 1;
    }

    /* Check kind is consistent with name */
    switch(func.info->kind) {
    case(XC_EXCHANGE):
      strcpy(kind, "_x_");
      break;

    case(XC_CORRELATION):
      strcpy(kind, "_c_");
      break;

    case(XC_EXCHANGE_CORRELATION):
      strcpy(kind, "_xc_");
      break;

    case(XC_KINETIC):
      strcpy(kind, "_k_");
      break;

    default:
      fprintf(stderr,"Kind %i not handled.\n",func.info->kind);
      return 2;
    }
    p = strstr(fname, kind);
    if(p == NULL)
      printf("Functional %i '%s' name may be inconsistent with its kind '%s'.\n",func_id, fname, kind);

    /* check if hybrid is initialized */
    if(strncmp(fname, "hyb_", 4) == 0){
      if(func.hyb_number_terms == 0)
        printf("Hybrid does not seem to be initialized\n");
    }

    /* Check family is consistent with name */
    family[0] = '\0';
    if(xc_hyb_type(&func) != XC_HYB_NONE)
      strcpy(family, "hyb_");

    switch(func.info->family) {
    case(XC_FAMILY_LDA):
      strcat(family, "lda_");
      break;

    case(XC_FAMILY_GGA):
      strcat(family, "gga_");
      break;

    case(XC_FAMILY_MGGA):
      strcat(family, "mgga_");
      break;

    default:
      fprintf(stderr, "Family %i not handled.\n", func.info->family);
      return 2;
    }

    p = strstr(fname, family);
    if(p != fname)
      printf("Functional %i '%s' name may be inconsistent with its family '%s'.\n",func_id, fname, family);

    /* Check non-local correlation parameters */
    {
      double delta, b, C;
      xc_hyb_vdw_vv10_coef(&func, &delta, &b, &C);

      if(func.hyb_type[0] == XC_HYB_VDW_VV10) {
        if(b == 0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero b.\n",func_id, fname);
        if(C == 0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero C.\n",func_id, fname);
      }
    }

    /* Free memory */
    free(fname);
    xc_func_end(&func);
  }

  /* Free memory */
  free(flist);

  return 0;
}

