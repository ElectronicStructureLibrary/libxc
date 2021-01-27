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

void print_hybrid(xc_func_type *func) {
  double alpha;
  assert(func->hyb_type[0] == XC_HYB_FOCK);
  alpha=func->hyb_coeff[0];
  printf("\nThis is a global hybrid functional with %4.1f%% of exact exchange.\n",alpha*100);
}

void print_rangesep(xc_func_type *func) {
  /* Determine amount of sr and lr exchange */
  double alpha, beta, omega;
  assert(func->hyb_type[0] == XC_HYB_ERF_SR || func->hyb_type[0] == XC_HYB_YUKAWA_SR);
  omega=func->hyb_omega[0];
  beta=func->hyb_coeff[0];
  if(func->hyb_number_terms>1) {
    assert(func->hyb_type[1] == XC_HYB_FOCK);
    alpha=func->hyb_coeff[1];
  } else {
    alpha=0.0;
  }
  printf("\nFunctional has a range-separation constant % .3f,\n",omega);
  printf("and %4.1f%% short-range and %4.1f%% long-range exact exchange,\n",(alpha+beta)*100,(alpha)*100);
}

void print_vdw(xc_func_type *func) {
  double Zab;
  assert(func->hyb_type[0] == XC_HYB_VDW_DF);
  Zab = func->hyb_coeff[0];
  printf("\nThis is a vdW-DF functional with Z_ab=%8.4f.\n", Zab);
}

int main(int argc, char **argv)
{
  int i, func_id, error, npar, hybrid_type;
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

  /* Print out hybrid exchange info */
  hybrid_type = xc_hyb_type(&func);
  switch(hybrid_type){
  case XC_HYB_SEMILOCAL:
    printf("\nThis is a semi-local functional with no exact exchange.\n");
    break;
  case XC_HYB_HYBRID:
    print_hybrid(&func);
    break;
  case XC_HYB_CAM:
    printf("\nThis is a range-separated functional, based on the error function kernel.\n");
    print_rangesep(&func);
    break;
  case XC_HYB_CAMY:
    printf("\nThis is a range-separated functional, based on the Yukawa kernel.\n");
    print_rangesep(&func);
    break;
  case XC_HYB_CAMG:
    printf("\nThis is a range-separated functional, based on the Gaussian attenuation scheme.\n");
    print_rangesep(&func);
    break;
  case XC_HYB_DOUBLE_HYBRID:
    printf("\nThis is a double hybrid functional.\n");
    break;
  case XC_HYB_VDW:
    printf("\nThis is a van der Waals functional.\n");
    print_vdw(&func);
    break;
  case XC_HYB_MIXTURE:
    printf("\nThis is a complicated functional ;)\n");
    break;

  default:
    printf("\nThis is functional of unknown type, please report this problem to the libxc tracker!\n");
    break;
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

  if(hybrid_type != XC_HYB_NONE) {
    printf("\nTo use the functional, the calling program must add in the following terms:\n");
    for(i=0; i<func.hyb_number_terms; i++){
      switch(func.hyb_type[i]){
      case XC_HYB_FOCK:
        printf("  *) Hartree-Fock exchange, weight = % .3f\n", func.hyb_coeff[i]);
        break;
      case XC_HYB_ERF_SR:
        printf("  *) Short-range (erf) exact exchange, weight = % .3f omega = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_YUKAWA_SR:
        printf("  *) Short-range (Yukawa) exact exchange, weight = % .3f omega = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_GAUSSIAN_SR:
        printf("  *) Short-range (Gau) exact exchange, weight = % .3f omega = %.3f\n", func.hyb_coeff[i], func.hyb_omega[i]);
        break;
      case XC_HYB_PT2:
        printf("  *) Second-order perturbation theory, weight = % .3f\n", func.hyb_coeff[i]);
        break;
      }
    }
  }

  /* Free memory */
  xc_func_end(&func);
  libxc_free(fname);

  return 0;
}
