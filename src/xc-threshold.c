/*
 Copyright (C) 2017 M.J.T. Oliveira

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <xc.h>

#define MIN_DENS 1.e-32

/*----------------------------------------------------------*/
void print_header(xc_func_type *func)
{

  printf("# Functional: %s\n#", func->info->name);

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

  printf("     zk    ");

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

  if (func->info->family == XC_FAMILY_LDA) {
    if (func->nspin == 1) {
      printf("    v3rho3 ");
    } else {
      printf("  v3rhoaaa   v3rhoaab   v3rhoabb   v3rhobbb ");
    }
  }

  printf("\n");
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  xc_func_type func, tf;
  int id, testcase, nspin;
  double a1, a2, b1, b2, delta;
  double r;
  double rho1a, rho1b, rho2a, rho2b;
  double grho1a, grho1b, grho2a, grho2b;
  double lapl1a, lapl1b, lapl2a, lapl2b;
  double tau1a, tau1b, tau2a, tau2b;
  double rho[2], sigma[3], lapl[2], tau[2];

  double zk;
  double vrho[2], vsigma[3], vlapl[2], vtau[2];
  double v2rho2[3], v2rhosigma[6], v2rholapl[3], v2rhotau[3], v2sigma2[6], v2sigmalapl[6],
      v2sigmatau[6], v2lapl2[3], v2lapltau[3], v2tau2[3];
  double v3rho3[4];

  /* Read input arguments and do some checks */
  if(argc != 3){
    printf("Usage:\n%s funct testcase\n", argv[0]);
    return 1;
  }
  id = atoi(argv[1]);
  testcase = atoi(argv[2]);

  if (testcase < 1 || testcase > 6) {
    fprintf(stderr, "Testcase %i not supported by xc-threshold.\nEnding program.\n", testcase);
    exit(1);
  }

  /* Initialize functional */
  nspin = testcase == 1 || testcase == 4 ? 1 : 2;
  if (xc_func_init(&func, id, nspin) != 0) {
    fprintf(stderr, "Functional '%d' not found\n", id);
    exit(1);
  }
  xc_func_set_dens_threshold(&func, MIN_DENS/10.0);


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
    rho[0] = rho1a + rho2a;
    sigma[0] = (grho1a + grho2a) * (grho1a + grho2a);
    lapl[0] = lapl1a + lapl2a;
    tau[0] = tau1a + tau2a;

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
    rho[1] = rho1b + rho2b;
    sigma[1] = (grho1a + grho2a) * (grho1b + grho2b);
    sigma[2] = (grho1b + grho2b) * (grho1b + grho2b);
    lapl[1] = lapl1b + lapl2b;
    tau[1] = tau1b + tau2b;

    /* Print densities and its derivatives */
    printf(" % 10.2e", rho[0]);
    if (nspin == 2) printf(" % 10.2e", rho[1]);
    if (func.info->family != XC_FAMILY_LDA) {
      printf(" % 10.2e", sigma[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e", sigma[1], sigma[2]);
    }
    if (func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA) {
      printf(" % 10.2e", lapl[0]);
      if (nspin == 2) printf(" % 10.2e", lapl[1]);
      printf(" % 10.2e", tau[0]);
      if (nspin == 2) printf(" % 10.2e", tau[1]);
    }

    /* Evaluate functional */
    switch(func.info->family) {
    case XC_FAMILY_LDA:
      xc_lda(&func, 1, rho, &zk, vrho, v2rho2, v3rho3);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga(&func, 1, rho, sigma,
	     &zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, NULL, NULL, NULL, NULL);
      break;

    case XC_FAMILY_MGGA:
    case XC_FAMILY_HYB_MGGA:
      xc_mgga(&func, 1, rho, sigma, lapl, tau,
	      &zk, vrho, vsigma, vlapl, vtau,
	      v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau, v2sigmalapl, v2sigmatau, v2lapltau);
      break;
    }

    /* Print results */
    printf(" % 10.2e", zk);

    printf(" % 10.2e", vrho[0]);
    if (nspin == 2) printf(" % 10.2e", vrho[1]);
    if (func.info->family != XC_FAMILY_LDA) {
      printf(" % 10.2e", vsigma[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e", vsigma[1], vsigma[2]);
    }
    if (func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA) {
      printf(" % 10.2e", vlapl[0]);
      if (nspin == 2) printf(" % 10.2e", vlapl[1]);
      printf(" % 10.2e", vtau[0]);
      if (nspin == 2) printf(" % 10.2e", vtau[1]);
    }

    printf(" % 10.2e", v2rho2[0]);
    if (nspin == 2) printf(" % 10.2e % 10.2e", v2rho2[1], v2rho2[2]);
    if (func.info->family != XC_FAMILY_LDA) {
      printf("  % 11.2e", v2sigma2[0]);
      if (nspin == 2) printf(" % 11.2e % 11.2e % 11.2e % 11.2e % 11.2e",
                             v2sigma2[1], v2sigma2[2], v2sigma2[3], v2sigma2[4], v2sigma2[5]);
    }
    if (func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA) {
      printf(" % 10.2e", v2lapl2[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e", v2lapl2[1], v2lapl2[2]);
      printf(" % 10.2e", v2tau2[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e", v2tau2[1], v2tau2[2]);
    }
    if (func.info->family != XC_FAMILY_LDA) {
      printf(" % 13.2e", v2rhosigma[0]);
      if (nspin == 2) printf(" % 13.2e % 13.2e % 13.2e % 13.2e % 13.2e",
                             v2rhosigma[1], v2rhosigma[2], v2rhosigma[3], v2rhosigma[4], v2rhosigma[5]);
    }
    if (func.info->family == XC_FAMILY_MGGA || func.info->family == XC_FAMILY_HYB_MGGA) {
      printf("  % 11.2e", v2rholapl[0]);
      if (nspin == 2) printf(" % 11.2e % 11.2e", v2rholapl[1], v2rhosigma[2]);
      printf("  % 10.2e", v2rhotau[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e", v2rhotau[1], v2rhotau[2]);
      printf(" % 14.2e", v2sigmalapl[0]);
      if (nspin == 2) printf(" % 14.2e % 14.2e % 14.2e % 14.2e % 14.2e",
                             v2sigmalapl[1], v2sigmalapl[2], v2sigmalapl[3], v2sigmalapl[4], v2sigmalapl[5]);
      printf("  % 13.2e", v2sigmatau[0]);
      if (nspin == 2) printf(" % 13.2e % 13.2e % 13.2e % 13.2e % 13.2e",
                             v2sigmatau[1], v2sigmatau[2], v2sigmatau[3], v2sigmatau[4], v2sigmatau[5]);
      printf("  % 11.2e", v2lapltau[0]);
      if (nspin == 2) printf(" % 11.2e % 11.2e", v2lapltau[1], v2lapltau[2]);
    }

    if (func.info->family == XC_FAMILY_LDA) {
      printf(" % 10.2e", v3rho3[0]);
      if (nspin == 2) printf(" % 10.2e % 10.2e % 10.2e", v3rho3[1], v3rho3[2], v3rho3[3]);
    }
    printf("\n");

    rho1a /= delta;
  } while (rho[0] > MIN_DENS && (rho[1] > MIN_DENS || rho[1] == 0.0));

  if (testcase > 3) xc_func_end(&tf);
  xc_func_end(&func);

  return 0;
}

