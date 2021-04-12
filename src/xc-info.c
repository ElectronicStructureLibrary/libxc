/*
 Copyright (C) 2014-2019 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <ctype.h>
#include "util.h"

#define is_mgga(id)   ((id) == XC_FAMILY_MGGA)
#define is_gga(id)    ((id) == XC_FAMILY_GGA || is_mgga(id))
#define is_lda(id)    ((id) == XC_FAMILY_LDA ||  is_gga(id))

int main(int argc, char **argv)
{
  int i, func_id, error, npar;
  xc_func_type func;
  char *fname;

  if(argc!=2) {
    printf("Usage: %s [ func_id | func_name ]\n",argv[0]);
    return 1;
  }

  /* Print libxc information */
  printf("libxc version %s\n",xc_version_string());
  printf("%s\n", xc_reference());
  printf("doi: %s\n", xc_reference_doi());
  printf("\n");

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

  printf("\nDefault thresholds:\n");
  printf("density: %e\n",func.dens_threshold);
  printf("   zeta: %e\n",func.zeta_threshold);
  if(is_gga(func.info->family))
    printf("  sigma: %e\n",func.sigma_threshold);
  if(is_mgga(func.info->family))
    printf("    tau: %e\n",func.tau_threshold);

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

  /* Print out hybrid exchange info */
  if(func.hyb_number_terms > 0){
    printf("\nTo use the functional, the calling program must add the following terms:\n");
    for(i=0; i<func.hyb_number_terms; i++){
      printf("  %d:", i);
      switch(func.hyb_type[i]){
      case XC_HYB_FOCK:
        printf("  Hartree-Fock exchange, weight = % .3f\n", func.hyb_params[i][0]);
        break;
      case XC_HYB_PT2:
        printf("  Second-order perturbation theory, weight = %.3f", func.hyb_params[i][0]);
        break;
      case XC_HYB_ERF_SR:
        printf("  Short-range (erf) exact exchange, weight = % .3f omega = %.3f\n",
               func.hyb_params[i][0], func.hyb_params[i][1]);
        break;
      case XC_HYB_YUKAWA_SR:
        printf("  Short-range (Yukawa) exact exchange, weight = % .3f omega = %.3f\n",
               func.hyb_params[i][0], func.hyb_params[i][1]);
        break;
      case XC_HYB_GAUSSIAN_SR:
        printf("  Short-range (Gaussian) exact exchange, weight = % .3f omega = %.3f\n",
               func.hyb_params[i][0], func.hyb_params[i][1]);
        break;
      case XC_HYB_VDW_DF:
        printf("  van der Waals (DF) with Zab = %.3f", func.hyb_params[i][0]);
        break;
      case XC_HYB_VDW_VV10:
        printf("  van der Waals (VV10) with b = %.3f and c = %.3f",
               func.hyb_params[i][0], func.hyb_params[i][1]);
        break;
      default:
        printf("\nThis term of unknown type, please report this problem to the libxc tracker!\n");
        break;
      }
    }
  }

  /* Free memory */
  xc_func_end(&func);
  libxc_free(fname);

  return 0;
}
