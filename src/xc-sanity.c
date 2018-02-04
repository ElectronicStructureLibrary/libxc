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
  char *fname, *p;

  /* Get list of available functionals */
  N=xc_number_of_functionals();
  flist=malloc(N*sizeof(int));
  xc_available_functional_numbers(flist);

  /* Loop over functionals */
  for(i=0;i<N;i++) {
    /* Functional id and name */
    int func_id;
    char *fname;
    /* Kind and family of functional */
    char kind[5], family[10];
    /* Hybrid parameters */
    double omega, alpha, beta;

    //printf("Checking functional %i -> %i\n",i,flist[i]);
    
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
      strcpy(kind,"_x_");
      break;
      
    case(XC_CORRELATION):
      strcpy(kind,"_c_");
      break;

    case(XC_EXCHANGE_CORRELATION):
      strcpy(kind,"_xc_");
      break;

    case(XC_KINETIC):
      strcpy(kind,"_k_");
      break;

    default:
      fprintf(stderr,"Kind %i not handled.\n",func.info->kind);
      return 2;
    }
    p=strstr(fname,kind);
    if(p == NULL)
      printf("Functional %i '%s' name may be inconsistent with its kind '%s'.\n",func_id, fname, kind);

    /* Check family is consistent with name */
    switch(func.info->family) {
    case(XC_FAMILY_LDA):
      strcpy(family,"lda_");
      break;
      
    case(XC_FAMILY_GGA):
      strcpy(family,"gga_");
      break;

    case(XC_FAMILY_MGGA):
      strcpy(family,"mgga_");
      break;

    case(XC_FAMILY_HYB_GGA):
      strcpy(family,"hyb_gga_");
      break;

    case(XC_FAMILY_HYB_MGGA):
      strcpy(family,"hyb_mgga_");
      break;

    default:
      fprintf(stderr,"Family %i not handled.\n",func.info->family);
      return 2;
    }
    p=strstr(fname,family);
    if(p != fname)
      printf("Functional %i '%s' name may be inconsistent with its family '%s'.\n",func_id, fname, family);
    
    /* Check hybrid parameters */
    xc_hyb_cam_coef(&func,&omega,&alpha,&beta);
    switch(func.info->family) {
    case(XC_FAMILY_LDA):
    case(XC_FAMILY_GGA):
    case(XC_FAMILY_MGGA):
      if(omega!=0.0)
        printf("Functional %i '%s' is supposed to be pure but has non-zero omega.\n",func_id, fname);
      if(alpha!=0.0)
        printf("Functional %i '%s' is supposed to be pure but has non-zero alpha.\n",func_id, fname);
      if(beta!=0.0)
        printf("Functional %i '%s' is supposed to be pure but has non-zero beta.\n",func_id, fname);
      break;

    case(XC_FAMILY_HYB_GGA):
    case(XC_FAMILY_HYB_MGGA):
      {
        /* Range separation? */
        int rangesep=0;
        if(func.info->flags & XC_FLAGS_HYB_CAM)
          rangesep++;
        if(func.info->flags & XC_FLAGS_HYB_CAMY)
          rangesep++;
        if(func.info->flags & XC_FLAGS_HYB_LC)
          rangesep++;
        if(func.info->flags & XC_FLAGS_HYB_LCY)
          rangesep++;

        switch(rangesep) {
        case(0):
          /* Conventional hybrid */
          if(alpha == 0.0)
            printf("Functional %i '%s' hybrid parameter hasn't been properly initialized.\n",func_id, fname);
          if(omega != 0.0 || beta != 0.0)
            printf("Functional %i '%s' erroneusly has range separated parameters.\n",func_id, fname);
          break;
          
        case(1):
          /* Range-separated hybrid */
          if(omega == 0.0)
            printf("Functional %i '%s' range separation constant hasn't been properly initialized.\n",func_id, fname);
          if(alpha == 0.0 && beta == 0.0)
            printf("Functional %i '%s' exact exchange parameters haven't been properly initialized.\n",func_id, fname);
          break;

        default:
          fprintf(stderr,"Functional %i '%s' has conflicting range separation flags set.\n",func_id,fname);
          return 3;
        }
      }
    }

    /* Check non-local correlation parameters */
    {
      double b, C;
      xc_nlc_coef(&func, &b, &C);
      
      if(func.info->flags & XC_FLAGS_VV10) {
        if(b==0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero b.\n",func_id, fname);
        if(C==0.0)
          printf("Functional %i '%s' is supposed to have VV10 but has zero b.\n",func_id, fname);
      } else {
        if(b != 0.0)
          printf("Functional %i '%s' isn't supposed to long-range correlation but has non-zero b.\n",func_id, fname);
        if(C != 0.0)
          printf("Functional %i '%s' isn't supposed to long-range correlation but has non-zero C.\n",func_id, fname);
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

