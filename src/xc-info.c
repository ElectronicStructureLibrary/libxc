/*
 Copyright (C) 2014-2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <ctype.h>
#include "util.h"

int main(int argc, char **argv)
{
  int i, func_id, error, npar, hybrid_type;
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
  printf("%10s: %-20i\n%10s: %-25s\n", "func_id", func_id, "name", fname);
  printf("%10s: %-20s\n%10s: %-25s\n", "family", get_family(&func), "kind", get_kind(&func));
  printf("%10s: %s\n", "comment", func.info->name);

  /* Print out hybrid exchange info */
  hybrid_type = xc_hyb_type(&func);
  switch(hybrid_type){
  case XC_HYB_SEMILOCAL:
    printf("\nThis is a semi-local functional.\n");
    break;
  case XC_HYB_HYBRID:
    printf("\nThis is a hybrid functional.\n");
    break;
  case XC_HYB_CAM:
    printf("\nThis is a CAM functional.\n");
    break;
  case XC_HYB_CAMY:
    printf("\nThis is a CAMY functional.\n");
    break;
  case XC_HYB_CAMG:
    printf("\nThis is a CAMG functional.\n");
    break;
  case XC_HYB_DOUBLE_HYBRID:
    printf("\nThis is a double hybrid functional.\n");
    break;
  case XC_HYB_MIXTURE:
    printf("\nThis is a complicated functional ;)\n");
    break;

  default:
    printf("\nThis is functional of unknown type, please report this problem to the libxc tracker!\n");
    break;
  }

  if(hybrid_type != XC_HYB_NONE) {
    printf("The calling program must add the following terms:\n");
    for(i=0; i<func.hyb_number_terms; i++){
      switch(func.hyb_type[i]){
      case XC_HYB_FOCK:
        printf("  *) Fock integral: mixing = %.3f\n", func.hyb_coeff[i]);
        break;
      case XC_HYB_ERF_SR:
        printf("  *) Short range Fock (erf): mixing = %.3f; mu = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_YUKAWA_SR:
        printf("  *) Short range Fock (Yukawa): mixing = %.3f; mu = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_GAUSSIAN_SR:
        printf("  *) Short range Fock (gaussian): mixing = %.3f; mu = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_PT2:
        printf("  *) Second order perturbation theory: mixing = %.3f", func.hyb_coeff[i]);
        break;

      }
    }
  }

  if(hybrid_type == XC_HYB_CAM || hybrid_type == XC_HYB_CAMY) {
    /* Determine amount of sr and lr exchange */
    double alpha, beta, omega;
    assert(func.hyb_type[0] == XC_HYB_ERF_SR || func.hyb_type[0] == XC_HYB_YUKAWA_SR);
    omega=func.hyb_omega[0];
    beta=func.hyb_coeff[0];
    assert(func.hyb_type[1] == XC_HYB_FOCK);
    alpha=func.hyb_coeff[1];
    printf("\nFunctional has a range-separation constant % .3f,\n",omega);
    printf("and %4.1f%% short-range and %4.1f%% long-range exact exchange,\n",(alpha+beta)*100,(alpha)*100);
  }

  printf("\nReference(s):\n");
  for(i = 0; i < 5; i++){
    if(func.info->refs[i] == NULL) break;
    printf("  *) %s\n", func.info->refs[i]->ref);
    if(strlen(func.info->refs[i]->doi) > 0){
       printf("     doi: %s\n", func.info->refs[i]->doi);
    }
  }

  printf("\nCompilation has support for:\n");
  if(func.info->flags & XC_FLAGS_HAVE_EXC)
    printf("  *) energy\n");
  if(func.info->flags & XC_FLAGS_HAVE_VXC)
    printf("  *) first derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_FXC)
    printf("  *) second derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_KXC)
    printf("  *) third derivative\n");
  if(func.info->flags & XC_FLAGS_HAVE_KXC)
    printf("  *) fourth derivative\n");

  printf("\nDefault density threshold: %e\n",func.dens_threshold);

  /* Query parameters */
  npar = xc_func_info_get_n_ext_params(func.info);
  if(npar > 0) {
    printf("\nFunctional has %i external parameters:\n",npar);
    printf("%3s %13s %8s %s\n","idx","value","name","description");
    for(i = 0; i < npar; i++)
      printf("%3i % e %8s %s\n", i,
             xc_func_info_get_ext_params_default_value(func.info, i),
             xc_func_info_get_ext_params_name(func.info, i),
             xc_func_info_get_ext_params_description(func.info, i));
  } else {
    printf("\nFunctional has no external parameters.\n");
  }

  /* Free memory */
  xc_func_end(&func);
  libxc_free(fname);

  return 0;
}
