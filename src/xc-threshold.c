/*
 Copyright (C) 2017 M.J.T. Oliveira

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <xc.h>
#include "util.h"

#define MIN_DENS 1.e-32

typedef struct {
  /* Input */
  double rho[2];         /* rhoa, rhob */
  double sigma[3];       /* sigmaaa, sigmaab, sigmabb */
  double lapl[2];        /* lapla, laplb */
  double tau[2];         /* taua, taub */

  /* Energy */
  double zk;             /* energy density per unit particle */

  /* First derivatives */
  double vrho[2];        /* vrhoa, vrhob */
  double vsigma[3];      /* vsigmaaa, vsigmaab, vsigmabb */
  double vlapl[2];       /* vlapla, vlaplb */
  double vtau[2];        /* vtaua, vtaub */

  /* Second derivatives */
  double v2rho2[3];      /* v2rhoa2, v2rhoab, v2rhob2 */
  double v2rhosigma[6];  /* v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb
	                          v2rhobsigmaaa, v2rhobsigmaab, v2rhobsigmabb */
  double v2rholapl[3];   /* */
  double v2rhotau[3];    /* */
  double v2sigma2[6];    /* v2sigmaaa2, v2sigmaaaab, v2sigmaaabb
			                      v2sigmaab2, v2sigmaabbb, v2sigmabb2 */
  double v2sigmalapl[6]; /* v2sigmaaalapla, v2sigmaaalaplb,
			                      v2sigmaablapla, v2sigmaablaplb,
			                      v2sigmabblapla, v2sigmabblaplb */
  double v2sigmatau[6];  /* v2sigmaaataua, v2sigmaaataub,
			                      v2sigmaabtaua, v2sigmaabtaub,
			                      v2sigmabbtaua, v2sigmabbtaub */
  double v2lapl2[3];     /* v2lapla2, v2laplab, v2laplb2 */
  double v2lapltau[3];   /* */
  double v2tau2[3];      /* v2taua2, v2tauab, v2taub2 */

  /* Third derivatives */
  double v3rho3[4];      /* v3rhoaaa, v3rhoaab, v3rhoabb, v3rhobbb */

} xc_values_type;

/*----------------------------------------------------------*/
void init_values(xc_values_type *xc_values)
{
  int i;

  xc_values->zk = 0;

  for(i=0; i<2; i++){
    xc_values->rho[i]   = 0;
    xc_values->lapl[i]  = 0;
    xc_values->tau[i]   = 0;
    xc_values->vrho[i]  = 0;
    xc_values->vlapl[i] = 0;
    xc_values->vtau[i]  = 0;
  }

  for(i=0; i<3; i++){
    xc_values->sigma[i]     = 0;
    xc_values->vsigma[i]    = 0;
    xc_values->v2rho2[i]    = 0;
    xc_values->v2lapl2[i]   = 0;
    xc_values->v2tau2[i]    = 0;
    xc_values->v2rholapl[i] = 0;
    xc_values->v2rhotau[i]  = 0;
    xc_values->v2lapltau[i] = 0;
  }

  for(i=0; i<4; i++){
    xc_values->v3rho3[i]  = 0;
  }

  for(i=0; i<6; i++){
    xc_values->v2rhosigma[i]  = 0;
    xc_values->v2sigma2[i]    = 0;
    xc_values->v2sigmalapl[i] = 0;
    xc_values->v2sigmatau[i]  = 0;
  }
}

/*----------------------------------------------------------*/
void print_header(xc_func_type *func)
{

  printf("# Functional: %s    Family: %s    Kind: %s\n", func->info->name, get_family(func), get_kind(func));

  printf("#");
  if (func->nspin == 1) {
    printf("    rho   ");
  } else {
    printf("    rhoa       rhob  ");
  }
  if (func->info->family != XC_FAMILY_LDA) {
    if (func->nspin == 1) {
      printf("    sigma  ");
    } else {
      printf("   sigmaaa    sigmaab    sigmabb ");
    }
  }
  if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
    if (func->nspin == 1) {
      printf("     lapl       tau   ");
    } else {
      printf("     lapla      laplb      taua       taub  ");
    }
  }

  if (func->info->flags & XC_FLAGS_HAVE_EXC)
    printf("     zk    ");

  if (func->info->flags & XC_FLAGS_HAVE_VXC) {
    if (func->nspin == 1) {
      printf("     vrho  ");
    } else {
      printf("    vrhoa      vrhob   ");
    }
    if (func->info->family != XC_FAMILY_LDA) {
      if (func->nspin == 1) {
        printf("    vsigma ");
      } else {
        printf("  vsigmaaa   vsigmaab   vsigmabb ");
      }
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      if (func->nspin == 1) {
        printf("    vlapl       vtau  ");
      } else {
        printf("   vlapla    vlaplb      vtaua      vtaub   ");
      }
    }
  }

  if (func->info->flags & XC_FLAGS_HAVE_FXC) {
    if (func->nspin == 1) {
      printf("    v2rho2 ");
    } else {
      printf("  v2rhoaa    v2rhoab    v2rhobb  ");
    }
    if (func->info->family != XC_FAMILY_LDA) {
      if (func->nspin == 1) {
        printf("     v2sigma2 ");
      } else {
        printf("  v2sigmaaaaa v2sigmaaaab v2sigmaaabb v2sigmaabab v2sigmaabbb v2sigmabbbb");
      }
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      if (func->nspin == 1) {
        printf("  v2lapl2     v2tau2  ");
      } else {
        printf("  v2laplaa   v2laplab   v2laplbb   v2tauaa    v2tauab    v2taubb  ");
      }
    }
    if (func->info->family != XC_FAMILY_LDA) {
      if (func->nspin == 1) {
        printf("    v2rhosigma ");
      } else {
        printf("  v2rhoasigmaaa v2rhoasigmaab v2rhoasigmabb v2rhobsigmaaa v2rhobsigmaab v2rhobsigmabb");
      }
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      if (func->nspin == 1) {
        printf("  v2rholapl    v2rhotau     v2sigmalapl     v2sigmatau   v2lapltau ");
      } else {
        printf(" v2rhoalapla v2rhoalaplb v2rhoblaplb ");
        printf(" v2rhoataua v2rhoataub v2rhobtaub");
        printf("   v2sigmaaalapla v2sigmaaalaplb v2sigmaablapla v2sigmaablaplb v2sigmabblapla v2sigmabblaplb");
        printf(" v2sigmaaataua v2sigmaaataub v2sigmaabtaua v2sigmaabtaub v2sigmabbtaua v2sigmabbtaub");
        printf(" v2laplataua v2laplataub v2laplbtaub");
      }
    }

  }

  if (func->info->flags & XC_FLAGS_HAVE_KXC) {
    if (func->info->family == XC_FAMILY_LDA) {
      if (func->nspin == 1) {
        printf("    v3rho3 ");
      } else {
        printf("  v3rhoaaa   v3rhoaab   v3rhoabb   v3rhobbb ");
      }
    }
  }

  printf("\n");
}

/*----------------------------------------------------------*/
void print_values(xc_values_type *values, xc_func_type *func){

  /* Print densities and its derivatives */
  printf(" % 10.2e", values->rho[0]);
  if (func->nspin == 2) printf(" % 10.2e", values->rho[1]);
  if (func->info->family != XC_FAMILY_LDA) {
    printf(" % 10.2e", values->sigma[0]);
    if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->sigma[1], values->sigma[2]);
  }
  if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
    printf(" % 10.2e", values->lapl[0]);
    if (func->nspin == 2) printf(" % 10.2e", values->lapl[1]);
    printf(" % 10.2e", values->tau[0]);
    if (func->nspin == 2) printf(" % 10.2e", values->tau[1]);
  }

  /* Print results */
  if (func->info->flags & XC_FLAGS_HAVE_EXC)
    printf(" % 10.2e", values->zk);

  if (func->info->flags & XC_FLAGS_HAVE_VXC) {
    printf(" % 10.2e", values->vrho[0]);
    if (func->nspin == 2) printf(" % 10.2e", values->vrho[1]);
    if (func->info->family != XC_FAMILY_LDA) {
      printf(" % 10.2e", values->vsigma[0]);
      if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->vsigma[1], values->vsigma[2]);
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      printf(" % 10.2e", values->vlapl[0]);
      if (func->nspin == 2) printf(" % 10.2e", values->vlapl[1]);
      printf(" % 10.2e", values->vtau[0]);
      if (func->nspin == 2) printf(" % 10.2e", values->vtau[1]);
    }
  }

  if (func->info->flags & XC_FLAGS_HAVE_FXC) {
    printf(" % 10.2e", values->v2rho2[0]);
    if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->v2rho2[1], values->v2rho2[2]);
    if (func->info->family != XC_FAMILY_LDA) {
      printf("  % 11.2e", values->v2sigma2[0]);
      if (func->nspin == 2) printf(" % 11.2e % 11.2e % 11.2e % 11.2e % 11.2e",
                                     values->v2sigma2[1], values->v2sigma2[2],
                                     values->v2sigma2[3], values->v2sigma2[4], values->v2sigma2[5]);
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      printf(" % 10.2e", values->v2lapl2[0]);
      if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->v2lapl2[1], values->v2lapl2[2]);
      printf(" % 10.2e", values->v2tau2[0]);
      if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->v2tau2[1], values->v2tau2[2]);
    }
    if (func->info->family != XC_FAMILY_LDA) {
      printf(" % 13.2e", values->v2rhosigma[0]);
      if (func->nspin == 2) printf(" % 13.2e % 13.2e % 13.2e % 13.2e % 13.2e",
                                     values->v2rhosigma[1], values->v2rhosigma[2],
                                     values->v2rhosigma[3], values->v2rhosigma[4], values->v2rhosigma[5]);
    }
    if (func->info->family == XC_FAMILY_MGGA || func->info->family == XC_FAMILY_HYB_MGGA) {
      printf("  % 11.2e", values->v2rholapl[0]);
      if (func->nspin == 2) printf(" % 11.2e % 11.2e", values->v2rholapl[1], values->v2rhosigma[2]);
      printf("  % 10.2e", values->v2rhotau[0]);
      if (func->nspin == 2) printf(" % 10.2e % 10.2e", values->v2rhotau[1], values->v2rhotau[2]);
      printf(" % 14.2e", values->v2sigmalapl[0]);
      if (func->nspin == 2) printf(" % 14.2e % 14.2e % 14.2e % 14.2e % 14.2e",
                                     values->v2sigmalapl[1], values->v2sigmalapl[2],
                                     values->v2sigmalapl[3], values->v2sigmalapl[4], values->v2sigmalapl[5]);
      printf("  % 13.2e", values->v2sigmatau[0]);
      if (func->nspin == 2) printf(" % 13.2e % 13.2e % 13.2e % 13.2e % 13.2e",
                                     values->v2sigmatau[1], values->v2sigmatau[2],
                                     values->v2sigmatau[3], values->v2sigmatau[4], values->v2sigmatau[5]);
      printf("  % 11.2e", values->v2lapltau[0]);
      if (func->nspin == 2) printf(" % 11.2e % 11.2e", values->v2lapltau[1], values->v2lapltau[2]);
    }
  }

  if (func->info->flags & XC_FLAGS_HAVE_KXC) {
    if (func->info->family == XC_FAMILY_LDA) {
      printf(" % 10.2e", values->v3rho3[0]);
      if (func->nspin == 2) printf(" % 10.2e % 10.2e % 10.2e", values->v3rho3[1], values->v3rho3[2], values->v3rho3[3]);
    }
  }

  printf("\n");
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  xc_func_type func, tf;
  xc_values_type values;
  int id, testcase, nspin;
  double a1, a2, b1, b2, delta;
  double r;
  double default_threshold;
  double rho1a, rho1b, rho2a, rho2b;
  double grho1a, grho1b, grho2a, grho2b;
  double lapl1a, lapl1b, lapl2a, lapl2b;
  double tau1a, tau1b, tau2a, tau2b;

  double *pzk          = NULL;
  double *pvrho        = NULL;
  double *pvsigma      = NULL;
  double *pvlapl       = NULL;
  double *pvtau        = NULL;
  double *pv2rho2      = NULL;
  double *pv2rhosigma  = NULL;
  double *pv2rholapl   = NULL;
  double *pv2rhotau    = NULL;
  double *pv2sigma2    = NULL;
  double *pv2sigmalapl = NULL;
  double *pv2sigmatau  = NULL;
  double *pv2lapl2     = NULL;
  double *pv2lapltau   = NULL;
  double *pv2tau2      = NULL;
  double *pv3rho3      = NULL;

  /* Read input arguments and do some checks */
  if(argc != 3){
    printf("Usage:\n%s funct testcase\n", argv[0]);
    return 1;
  }
  id = atoi(argv[1]);
  testcase = atoi(argv[2]);

  switch (testcase) {
    case 1:
      printf("# Test case: %s", "atomic densities, spin-unpolarized");
      break;
    case 2:
      printf("# Test case: %s", "atomic densities, spin-polarized");
      break;
    case 3:
      printf("# Test case: %s", "atomic densities, spin-polarized, one empty spin-channel");
      break;
    case 4:
      printf("# Test case: %s", "HEG, spin-unpolarized");
      break;
    case 5:
      printf("# Test case: %s", "HEG, spin-polarized");
      break;
    case 6:
      printf("# Test case: %s", "HEG, spin-polarized, one empty spin-channel");
      break;
    default:
      fprintf(stderr, "Testcase %i not supported by xc-threshold.\nEnding program.\n", testcase);
      exit(1);
  }
  printf("\n");

  /* Initialize values */
  init_values(&values);

  /* Initialize functional */
  nspin = testcase == 1 || testcase == 4 ? 1 : 2;
  if (xc_func_init(&func, id, nspin) != 0) {
    fprintf(stderr, "Functional '%d' not found\n", id);
    exit(1);
  }
  default_threshold = func.dens_threshold;
  xc_func_set_dens_threshold(&func, MIN_DENS/10.0);
  if(func.info->flags & XC_FLAGS_HAVE_EXC){
    pzk = &values.zk;
  }
  if(func.info->flags & XC_FLAGS_HAVE_VXC){
    pvrho   = values.vrho;
    pvsigma = values.vsigma;
    pvlapl  = values.vlapl;
    pvtau   = values.vtau;
  }
  if(func.info->flags & XC_FLAGS_HAVE_FXC){
    pv2rho2      = values.v2rho2;
    pv2rhosigma  = values.v2rhosigma;
    pv2rholapl   = values.v2rholapl;
    pv2rhotau    = values.v2rhotau;
    pv2sigma2    = values.v2sigma2;
    pv2sigmalapl = values.v2sigmalapl;
    pv2sigmatau  = values.v2sigmatau;
    pv2lapl2     = values.v2lapl2;
    pv2lapltau   = values.v2lapltau;
    pv2tau2      = values.v2tau2;
  }
  if(func.info->flags & XC_FLAGS_HAVE_KXC){
    pv3rho3 = values.v3rho3;
  }

  /* The trial densities are build from two states (four in case of a spin-polarized system).
   * This avoids problems with MGGA's, as some terms vanish for one state systems.
   *
   * The following parameters control the exponential decay of the states densities:
   *  rho_1^a(r) = e^{- 2 a1 r}
   *  rho_1^b(r) = e^{- 2 b1 r}
   *  rho_2^a(r) = e^{- 2 a2 r}
   *  rho_2^b(r) = e^{- 2 b2 r}
   *
   * The sampling is determined by the delta parameter:
   *  r_{i+1} = r_i + log(delta)/a1
   */
  a1 = 1.0;  a2 = 0.75; b1 = 0.5; b2 = 0.25;
  delta = 8.0;

  /* Print header */
  print_header(&func);

  /* Initialize some values to zero so that we do not need to do it during the loop */
  rho1b = rho2b = 0.0;
  grho1a = grho2a = grho1b = grho2b = 0.0;
  lapl1a = lapl2a = lapl1b = lapl2b = 0.0;
  tau1a = tau2a = tau1b = tau2b = 0.0;

  /* Auxiliary TF functional */
  if (testcase > 3) xc_func_init(&tf, XC_LDA_K_TF, nspin);

  rho1a = 0.5; /* Initial value for state 1a density */
  do {
    r = -log(rho1a) / (2.0 * a1);

    /* Compute densities and derivatives.
     *
     * There are 6 test cases implemented:
     *  - case 1: Atomic densities, spin-unpolarized.
     *  - case 2: Atomic densities, spin-polarized.
     *  - case 3: Atomic densities, spin-polarized, one empty spin-channel.
     *  - case 4: HEG, spin-unpolarized.
     *  - case 5: HEG, spin-polarized.
     *  - case 6: HEG, spin-polarized, one empty spin-channel.
     */
    rho2a = exp(-2.0 * a2 * r);
    if (testcase == 1 || testcase == 2 || testcase == 3) {
      grho1a = -2.0 * a1 * rho1a;
      grho2a = -2.0 * a2 * rho2a;
      lapl1a = 4.0 * a1 * (a1 + 1.0 / r) * rho1a;
      lapl2a = 4.0 * a2 * (a2 + 1.0 / r) * rho2a;
      tau1a = 4.0 * a1 * a1 * rho1a;
      tau2a = 4.0 * a2 * a2 * rho2a;
    }
    if (testcase > 3) {
      xc_lda_exc(&tf, 1, &rho1a, &tau1a);
      xc_lda_exc(&tf, 1, &rho2a, &tau2a);
    }
    if (testcase == 2 || testcase == 5) {
      rho1b = exp(-2.0 * b1 * r);
      rho2b = exp(-2.0 * b2 * r);
    }
    if (testcase == 2) {
      grho1b = -2.0 * b1 * rho1b;
      grho2b = -2.0 * b2 * rho2b;
      lapl1b = 4.0 * b1 * (b1 + 1.0 / r) * rho1b;
      lapl2b = 4.0 * b2 * (b2 + 1.0 / r) * rho2b;
      tau1b = 4.0 * b1 * b1 * rho1b;
      tau2b = 4.0 * b2 * b2 * rho2b;
    }
    if (testcase == 5) {
      xc_lda_exc(&tf, 1, &rho1b, &tau1b);
      xc_lda_exc(&tf, 1, &rho2b, &tau2b);
    }

    values.rho[0]   = rho1a + rho2a;
    values.rho[1]   = rho1b + rho2b;
    values.sigma[0] = (grho1a + grho2a) * (grho1a + grho2a);
    values.sigma[1] = (grho1a + grho2a) * (grho1b + grho2b);
    values.sigma[2] = (grho1b + grho2b) * (grho1b + grho2b);
    values.lapl[0]  = lapl1a + lapl2a;
    values.lapl[1]  = lapl1b + lapl2b;
    values.tau[0]   = tau1a + tau2a;
    values.tau[1]   = tau1b + tau2b;

    if (default_threshold > values.rho[0]) {
      printf("# default threshold\n");
      default_threshold = MIN_DENS/10.;
    }

    /* Evaluate functional */
    switch(func.info->family) {
    case XC_FAMILY_LDA:
      xc_lda(&func, 1, values.rho, pzk, pvrho, pv2rho2, pv3rho3);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga(&func, 1, values.rho, values.sigma,
	     pzk, pvrho, pvsigma, pv2rho2, pv2rhosigma, pv2sigma2, NULL, NULL, NULL, NULL);
      break;

    case XC_FAMILY_MGGA:
    case XC_FAMILY_HYB_MGGA:
      xc_mgga(&func, 1, values.rho, values.sigma, values.lapl, values.tau,
	      pzk, pvrho, pvsigma, pvlapl, pvtau,
	      pv2rho2, pv2sigma2, pv2lapl2, pv2tau2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigmalapl, pv2sigmatau, pv2lapltau);
      break;
    }

    print_values(&values, &func);

    rho1a /= delta;
  } while (values.rho[0] > MIN_DENS && (values.rho[1] > MIN_DENS || values.rho[1] == 0.0));

  if (testcase > 3) xc_func_end(&tf);
  xc_func_end(&func);

  return 0;
}

