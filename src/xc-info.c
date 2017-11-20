/*
 Copyright (C) 2014 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <ctype.h>
#include "util.h"

int main(int argc, char **argv) {
  int i, func_id, error;
  xc_func_type func;
  char *fname;

  if(argc!=2) {
    printf("Usage: %s [ func_id | func_name ]\n",argv[0]);
    return 1;
  }

  /* Is functional defined by a string constant? */
  if(isalpha(argv[1][0]))
    func_id = xc_functional_get_number(argv[1]);
  else
    func_id = atoi(argv[1]);

  /* Initialize functional */
  error = xc_func_init(&func, func_id, XC_UNPOLARIZED);
  if(error) {
    printf("Functional '%s' not found.\n", argv[1]);
    return 1;
  }

  /* Get functional name */
  fname = xc_functional_get_name(func_id);

  /* Print out info */
  printf("%10s: %-20i\t%10s: %-25s\n","func_id", func_id, "name", fname);
  printf("%10s: %-20s\t%10s: %-25s\n","family", get_family(&func), "kind", get_kind(&func));
  printf("%10s: %s\n","comment", func.info->name);

  /* Print out hybrid exchange info */
  if(func.info->family==XC_FAMILY_HYB_GGA || func.info->family==XC_FAMILY_HYB_MGGA) {
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

    if(rangesep) {
      double alpha, beta, omega;
      xc_hyb_cam_coef(&func,&omega,&alpha,&beta);
      printf("\nThis is a range-separated hybrid functional with range-separation constant % .3f,\n",omega);
      printf("and %4.1f%% short-range and %4.1f%% long-range exact exchange.\n",(alpha+beta)*100,(alpha)*100);
    } else {
      double alpha=xc_hyb_exx_coef(&func);
      printf("\nThis is a global hybrid functional with %4.1f%% of exact exchange.\n",alpha*100);
    }
  } else {
    if(func.info->kind == XC_EXCHANGE || func.info->kind == XC_EXCHANGE_CORRELATION)
      printf("\nThis is a pure functional with no exact exchange.\n");
  }
  
  printf("\nReference(s):\n");
  for(i=0; i<5; i++){
    if(func.info->refs[i]==NULL) break;
    printf("%s", func.info->refs[i]->ref);
    if(strlen(func.info->refs[i]->doi) > 0){
       printf(" (%s)", func.info->refs[i]->doi);
    }
    printf("\n");
  }

  /* Free memory */
  xc_func_end(&func);
  free(fname);

  return 0;
}
