/*
 Copyright (C) 2017 M.J.T. Oliveira
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include <xc.h>

/* Get debug info by uncommenting the following line */
/* #define XC_DEBUG */

typedef struct {
  /* Input */
  double rho[2];   /* rhoa, rhob */
  double sigma[3]; /* sigmaaa, sigmaab, sigmabb */
  double lapl[2];  /* lapla, laplb */
  double tau[2];   /* taua, taub */

  /* Output of xc-dimensions */
  double zk;
  double vrho[2];
  double vsigma[3];
  double vlapl[2];
  double vtau[2];
  double v2rho2[3];
  double v2rhosigma[6];
  double v2sigma2[6];
  double v2rholapl[4];
  double v2sigmalapl[6];
  double v2lapl2[3];
  double v2rhotau[4];
  double v2sigmatau[6];
  double v2lapltau[4];
  double v2tau2[3];
  double v3rho3[4];
  double v3rho2sigma[9];
  double v3rhosigma2[12];
  double v3sigma3[10];
  double v3rho2lapl[6];
  double v3rhosigmalapl[12];
  double v3sigma2lapl[12];
  double v3rholapl2[6];
  double v3sigmalapl2[9];
  double v3lapl3[4];
  double v3rho2tau[6];
  double v3rhosigmatau[12];
  double v3sigma2tau[12];
  double v3rholapltau[8];
  double v3sigmalapltau[12];
  double v3lapl2tau[6];
  double v3rhotau2[6];
  double v3sigmatau2[9];
  double v3lapltau2[6];
  double v3tau3[4];
  double v4rho4[5];
  double v4rho3sigma[12];
  double v4rho2sigma2[18];
  double v4rhosigma3[20];
  double v4sigma4[15];
  double v4rho3lapl[8];
  double v4rho2sigmalapl[18];
  double v4rhosigma2lapl[36];
  double v4sigma3lapl[20];
  double v4rho2lapl2[9];
  double v4rhosigmalapl2[18];
  double v4sigma2lapl2[18];
  double v4rholapl3[8];
  double v4sigmalapl3[12];
  double v4lapl4[5];
  double v4rho3tau[8];
  double v4rho2sigmatau[18];
  double v4rhosigma2tau[36];
  double v4sigma3tau[30];
  double v4rho2lapltau[12];
  double v4rhosigmalapltau[24];
  double v4sigma2lapltau[24];
  double v4rholapl2tau[12];
  double v4sigmalapl2tau[18];
  double v4lapl3tau[8];
  double v4rho2tau2[9];
  double v4rhosigmatau2[36];
  double v4sigma2tau2[18];
  double v4rholapltau2[12];
  double v4sigmalapltau2[18];
  double v4lapl2tau2[9];
  double v4rhotau3[8];
  double v4sigmatau3[12];
  double v4lapltau3[8];
  double v4tau4[5];
} xc_values_type;

void init_values(xc_values_type *xc_values) {
  memset(xc_values, 0, sizeof(xc_values_type));
}

/* Implemented test cases */
typedef enum {
  /* Atomic densities, spin-unpolarized */
  ATOMIC_UPOL,
  /* Atomic densities, spin-polarized */
  ATOMIC_POL,
  /* Atomic densities, spin-polarized with empty spin channel */
  ATOMIC_EMPTY,
  /* Homogeneous electron gas, spin-unpolarized */
  HEG_UPOL,
  /* Homogeneous electron gas, spin-polarized */
  HEG_POL,
  /* Homogeneous electron gas, spin-polarized with empty spin channel */
  HEG_EMPTY,
  /* Fallthrough */
  NOSUCH_TEST
} testcase_t;

void compute_input(xc_values_type *values, int *nspin, testcase_t testcase,
                   double r) {
  /* The trial densities are built from two states (four in case of a
   * spin-polarized system).  This avoids problems with MGGA's, as
   * some terms vanish for one state systems.
   *
   * The following parameters control the exponential decay of the states
   * densities: rho_1^a(r) = e^{- 2 a1 r} rho_1^b(r) = e^{- 2 b1 r} rho_2^a(r) =
   * e^{- 2 a2 r} rho_2^b(r) = e^{- 2 b2 r}
   */
  const double a1 = 1.0;
  const double a2 = 0.75;
  const double b1 = 0.5;
  const double b2 = 0.25;

  /* Input variables */
  double rho1a, rho1b, rho2a, rho2b;
  double grho1a, grho1b, grho2a, grho2b;
  double lapl1a, lapl1b, lapl2a, lapl2b;
  double tau1a, tau1b, tau2a, tau2b;

  /* Helper function */
  xc_func_type tf;

  /* Default values for computed stuff  */
  rho1b = rho2b = 0.0;
  grho1a = grho2a = grho1b = grho2b = 0.0;
  lapl1a = lapl2a = lapl1b = lapl2b = 0.0;
  tau1a = tau2a = tau1b = tau2b = 0.0;

  /* Compute atomic inputs. */
  rho1a = exp(-2.0 * a1 * r);
  rho2a = exp(-2.0 * a2 * r);
  rho1b = exp(-2.0 * b1 * r);
  rho2b = exp(-2.0 * b2 * r);

  grho1a = -2.0 * a1 * rho1a;
  grho2a = -2.0 * a2 * rho2a;
  grho1b = -2.0 * b1 * rho1b;
  grho2b = -2.0 * b2 * rho2b;

  lapl1a = 4.0 * a1 * (a1 + 1.0 / r) * rho1a;
  lapl2a = 4.0 * a2 * (a2 + 1.0 / r) * rho2a;
  lapl1b = 4.0 * b1 * (b1 + 1.0 / r) * rho1b;
  lapl2b = 4.0 * b2 * (b2 + 1.0 / r) * rho2b;

  tau1a = 4.0 * a1 * a1 * rho1a;
  tau2a = 4.0 * a2 * a2 * rho2a;
  tau1b = 4.0 * b1 * b1 * rho1b;
  tau2b = 4.0 * b2 * b2 * rho2b;

  /* If HEG, overwrite with HEG values */
  if (testcase > ATOMIC_EMPTY) {
    xc_func_init(&tf, XC_LDA_K_TF, XC_UNPOLARIZED);
    xc_lda_exc(&tf, 1, &rho1a, &tau1a);
    xc_lda_exc(&tf, 1, &rho2a, &tau2a);
    xc_lda_exc(&tf, 1, &rho1b, &tau1b);
    xc_lda_exc(&tf, 1, &rho2b, &tau2b);
  }

  /* Spin-restrict or zero */
  if (testcase == ATOMIC_EMPTY || testcase == HEG_EMPTY) {
    rho1b = 0.0;
    rho2b = 0.0;
    grho1b = 0.0;
    grho2b = 0.0;
    lapl1b = 0.0;
    lapl2b = 0.0;
    tau1b = 0.0;
    tau2b = 0.0;
  } else if (testcase == ATOMIC_UPOL || testcase == HEG_UPOL) {
    rho1b = rho1a;
    rho2b = rho2a;
    grho1b = grho1a;
    grho2b = grho2a;
    lapl1b = lapl1a;
    lapl2b = lapl1b;
  }

  /* Initialize libxc input and output values */
  init_values(values);
  /* Compute input */
  values->rho[0] = rho1a + rho2a;
  values->rho[1] = rho1b + rho2b;
  values->sigma[0] = (grho1a + grho2a) * (grho1a + grho2a);
  values->sigma[1] = (grho1a + grho2a) * (grho1b + grho2b);
  values->sigma[2] = (grho1b + grho2b) * (grho1b + grho2b);
  values->lapl[0] = lapl1a + lapl2a;
  values->lapl[1] = lapl1b + lapl2b;
  values->tau[0] = tau1a + tau2a;
  values->tau[1] = tau1b + tau2b;

  /* End helper */
  if (testcase > ATOMIC_EMPTY)
    xc_func_end(&tf);

  /* Store spin */
  *nspin = (testcase == ATOMIC_UPOL || testcase == HEG_UPOL) ? XC_UNPOLARIZED
                                                             : XC_POLARIZED;
}

size_t check_xc(int id, int nspin, xc_values_type values, double threshold) {
  size_t nfail = 0;
  xc_func_type func;

  double *prho = values.rho;
  double *psigma = values.sigma;
  double *plapl = values.lapl;
  double *ptau = values.tau;

  double *pzk = NULL;

  /* This is from gen-xc-pointers.sh */
  double *pvrho = NULL;
  double *pvsigma = NULL;
  double *pvlapl = NULL;
  double *pvtau = NULL;
  double *pv2rho2 = NULL;
  double *pv2rhosigma = NULL;
  double *pv2sigma2 = NULL;
  double *pv2rholapl = NULL;
  double *pv2sigmalapl = NULL;
  double *pv2lapl2 = NULL;
  double *pv2rhotau = NULL;
  double *pv2sigmatau = NULL;
  double *pv2lapltau = NULL;
  double *pv2tau2 = NULL;
  double *pv3rho3 = NULL;
  double *pv3rho2sigma = NULL;
  double *pv3rhosigma2 = NULL;
  double *pv3sigma3 = NULL;
  double *pv3rho2lapl = NULL;
  double *pv3rhosigmalapl = NULL;
  double *pv3sigma2lapl = NULL;
  double *pv3rholapl2 = NULL;
  double *pv3sigmalapl2 = NULL;
  double *pv3lapl3 = NULL;
  double *pv3rho2tau = NULL;
  double *pv3rhosigmatau = NULL;
  double *pv3sigma2tau = NULL;
  double *pv3rholapltau = NULL;
  double *pv3sigmalapltau = NULL;
  double *pv3lapl2tau = NULL;
  double *pv3rhotau2 = NULL;
  double *pv3sigmatau2 = NULL;
  double *pv3lapltau2 = NULL;
  double *pv3tau3 = NULL;
  double *pv4rho4 = NULL;
  double *pv4rho3sigma = NULL;
  double *pv4rho2sigma2 = NULL;
  double *pv4rhosigma3 = NULL;
  double *pv4sigma4 = NULL;
  double *pv4rho3lapl = NULL;
  double *pv4rho2sigmalapl = NULL;
  double *pv4rhosigma2lapl = NULL;
  double *pv4sigma3lapl = NULL;
  double *pv4rho2lapl2 = NULL;
  double *pv4rhosigmalapl2 = NULL;
  double *pv4sigma2lapl2 = NULL;
  double *pv4rholapl3 = NULL;
  double *pv4sigmalapl3 = NULL;
  double *pv4lapl4 = NULL;
  double *pv4rho3tau = NULL;
  double *pv4rho2sigmatau = NULL;
  double *pv4rhosigma2tau = NULL;
  double *pv4sigma3tau = NULL;
  double *pv4rho2lapltau = NULL;
  double *pv4rhosigmalapltau = NULL;
  double *pv4sigma2lapltau = NULL;
  double *pv4rholapl2tau = NULL;
  double *pv4sigmalapl2tau = NULL;
  double *pv4lapl3tau = NULL;
  double *pv4rho2tau2 = NULL;
  double *pv4rhosigmatau2 = NULL;
  double *pv4sigma2tau2 = NULL;
  double *pv4rholapltau2 = NULL;
  double *pv4sigmalapltau2 = NULL;
  double *pv4lapl2tau2 = NULL;
  double *pv4rhotau3 = NULL;
  double *pv4sigmatau3 = NULL;
  double *pv4lapltau3 = NULL;
  double *pv4tau4 = NULL;

  /* Initialize functional */
  if (xc_func_init(&func, id, nspin) != 0) {
    fprintf(stderr, "Error initializing functional id '%d'\n", id);
    exit(1);
  }

  /* Initialize pointers, this is from gen-xc-pointers.sh */
  if (func.info->flags & XC_FLAGS_HAVE_VXC) {
    pvrho = values.vrho;
    pvsigma = values.vsigma;
    pvlapl = values.vlapl;
    pvtau = values.vtau;
  }
  if (func.info->flags & XC_FLAGS_HAVE_FXC) {
    pv2rho2 = values.v2rho2;
    pv2rhosigma = values.v2rhosigma;
    pv2sigma2 = values.v2sigma2;
    pv2rholapl = values.v2rholapl;
    pv2sigmalapl = values.v2sigmalapl;
    pv2lapl2 = values.v2lapl2;
    pv2rhotau = values.v2rhotau;
    pv2sigmatau = values.v2sigmatau;
    pv2lapltau = values.v2lapltau;
    pv2tau2 = values.v2tau2;
  }
  if (func.info->flags & XC_FLAGS_HAVE_KXC) {
    pv3rho3 = values.v3rho3;
    pv3rho2sigma = values.v3rho2sigma;
    pv3rhosigma2 = values.v3rhosigma2;
    pv3sigma3 = values.v3sigma3;
    pv3rho2lapl = values.v3rho2lapl;
    pv3rhosigmalapl = values.v3rhosigmalapl;
    pv3sigma2lapl = values.v3sigma2lapl;
    pv3rholapl2 = values.v3rholapl2;
    pv3sigmalapl2 = values.v3sigmalapl2;
    pv3lapl3 = values.v3lapl3;
    pv3rho2tau = values.v3rho2tau;
    pv3rhosigmatau = values.v3rhosigmatau;
    pv3sigma2tau = values.v3sigma2tau;
    pv3rholapltau = values.v3rholapltau;
    pv3sigmalapltau = values.v3sigmalapltau;
    pv3lapl2tau = values.v3lapl2tau;
    pv3rhotau2 = values.v3rhotau2;
    pv3sigmatau2 = values.v3sigmatau2;
    pv3lapltau2 = values.v3lapltau2;
    pv3tau3 = values.v3tau3;
  }
  if (func.info->flags & XC_FLAGS_HAVE_LXC) {
    pv4rho4 = values.v4rho4;
    pv4rho3sigma = values.v4rho3sigma;
    pv4rho2sigma2 = values.v4rho2sigma2;
    pv4rhosigma3 = values.v4rhosigma3;
    pv4sigma4 = values.v4sigma4;
    pv4rho3lapl = values.v4rho3lapl;
    pv4rho2sigmalapl = values.v4rho2sigmalapl;
    pv4rhosigma2lapl = values.v4rhosigma2lapl;
    pv4sigma3lapl = values.v4sigma3lapl;
    pv4rho2lapl2 = values.v4rho2lapl2;
    pv4rhosigmalapl2 = values.v4rhosigmalapl2;
    pv4sigma2lapl2 = values.v4sigma2lapl2;
    pv4rholapl3 = values.v4rholapl3;
    pv4sigmalapl3 = values.v4sigmalapl3;
    pv4lapl4 = values.v4lapl4;
    pv4rho3tau = values.v4rho3tau;
    pv4rho2sigmatau = values.v4rho2sigmatau;
    pv4rhosigma2tau = values.v4rhosigma2tau;
    pv4sigma3tau = values.v4sigma3tau;
    pv4rho2lapltau = values.v4rho2lapltau;
    pv4rhosigmalapltau = values.v4rhosigmalapltau;
    pv4sigma2lapltau = values.v4sigma2lapltau;
    pv4rholapl2tau = values.v4rholapl2tau;
    pv4sigmalapl2tau = values.v4sigmalapl2tau;
    pv4lapl3tau = values.v4lapl3tau;
    pv4rho2tau2 = values.v4rho2tau2;
    pv4rhosigmatau2 = values.v4rhosigmatau2;
    pv4sigma2tau2 = values.v4sigma2tau2;
    pv4rholapltau2 = values.v4rholapltau2;
    pv4sigmalapltau2 = values.v4sigmalapltau2;
    pv4lapl2tau2 = values.v4lapl2tau2;
    pv4rhotau3 = values.v4rhotau3;
    pv4sigmatau3 = values.v4sigmatau3;
    pv4lapltau3 = values.v4lapltau3;
    pv4tau4 = values.v4tau4;
  }

  /* Set the threshold */
  xc_func_set_dens_threshold(&func, threshold);
  /* Evaluate functional */
  switch (func.info->family) {
  case XC_FAMILY_LDA:
    xc_lda(&func, 1, prho, pzk, pvrho, pv2rho2, pv3rho3, pv4rho4);
    break;
  case XC_FAMILY_GGA:
    xc_gga(&func, 1, prho, psigma, pzk, pvrho, pvsigma, pv2rho2, pv2rhosigma,
           pv2sigma2, pv3rho3, pv3rho2sigma, pv3rhosigma2, pv3sigma3, pv4rho4,
           pv4rho3sigma, pv4rho2sigma2, pv4rhosigma3, pv4sigma4);
    break;
  case XC_FAMILY_MGGA:
    xc_mgga(&func, 1, prho, psigma, plapl, ptau, pzk, pvrho, pvsigma, pvlapl,
            pvtau, pv2rho2, pv2rhosigma, pv2rholapl, pv2rhotau, pv2sigma2,
            pv2sigmalapl, pv2sigmatau, pv2lapl2, pv2lapltau, pv2tau2, pv3rho3,
            pv3rho2sigma, pv3rho2lapl, pv3rho2tau, pv3rhosigma2,
            pv3rhosigmalapl, pv3rhosigmatau, pv3rholapl2, pv3rholapltau,
            pv3rhotau2, pv3sigma3, pv3sigma2lapl, pv3sigma2tau, pv3sigmalapl2,
            pv3sigmalapltau, pv3sigmatau2, pv3lapl3, pv3lapl2tau, pv3lapltau2,
            pv3tau3, pv4rho4, pv4rho3sigma, pv4rho3lapl, pv4rho3tau,
            pv4rho2sigma2, pv4rho2sigmalapl, pv4rho2sigmatau, pv4rho2lapl2,
            pv4rho2lapltau, pv4rho2tau2, pv4rhosigma3, pv4rhosigma2lapl,
            pv4rhosigma2tau, pv4rhosigmalapl2, pv4rhosigmalapltau,
            pv4rhosigmatau2, pv4rholapl3, pv4rholapl2tau, pv4rholapltau2,
            pv4rhotau3, pv4sigma4, pv4sigma3lapl, pv4sigma3tau, pv4sigma2lapl2,
            pv4sigma2lapltau, pv4sigma2tau2, pv4sigmalapl3, pv4sigmalapl2tau,
            pv4sigmalapltau2, pv4sigmatau3, pv4lapl4, pv4lapl3tau, pv4lapl2tau2,
            pv4lapltau3, pv4tau4);
    break;

  default:
    fprintf(stderr, "Unsupported functional family %i!\n", func.info->family);
  }

  /* Check the output for NaNs */
  double *ptr = (double *)&values;
  for (int i = 0; i < sizeof(xc_values_type) / sizeof(double); i++)
    /* We have encountered an infinity or NaN */
    if (!isfinite(ptr[i])) {
#ifdef XC_DEBUG
      printf("values[%i] = %e\n", i, ptr[i]);
#endif
      nfail++;
    }

  /* Deallocate functional */
  xc_func_end(&func);

  return nfail;
}

double get_threshold(int id) {
  double thr;
  xc_func_type func;
  if (xc_func_init(&func, id, XC_POLARIZED) != 0) {
    fprintf(stderr, "Error initializing functional id '%d'\n", id);
    exit(1);
  }
  thr = func.dens_threshold;
  xc_func_end(&func);
  return thr;
}

int main(int argc, char *argv[]) {
  /* xc values */
  xc_values_type values;
  /* functional id */
  int id;
  /* spin type */
  int nspin;

  /* test case */
  int testcase;
  /* Current threshold */
  double threshold;

  /* logarithmic grid spacing */
  const double h = 0.025;
  /* density at end point */
  const double nmin = 1e-60;
  /* since we use an exponential density with zeta = 1.0, the
     necessary maximal radius is */
  const double rmax = -log(nmin) / 2;
  /* this leads to the number of radial points as */
  const double nrad = ceil(log(rmax) / h);

  /* radial point id */
  int irad;
  /* radius */
  double r;
  /* functional name */
  char *fname;

  /* Read in function id */
  if (argc != 2) {
    printf("Usage:\n%s funct\n", argv[0]);
    return 1;
  }

#if 0
  printf("logarithmic spacing h = %e yields %lu radial points\n", h,
         (long unsigned)nrad);
  printf("rmax = %e, nmin = %e\n", exp(h * (nrad - 1)),
         exp(-2 * exp(h * (nrad - 1))));
#endif
  
  /* Is functional defined by a string constant? */
  if (isalpha(argv[1][0]))
    id = xc_functional_get_number(argv[1]);
  else
    id = atoi(argv[1]);

  /* Get functional name */
  fname = xc_functional_get_name(id);
  printf("Functional id %i : %s\n", id, fname);

  /* Initial value for threshold */
  threshold = 1.0;
  while (1) {
    /* Number of failed data points */
    size_t nfail = 0;
    /* Run all test cases for the given threshold */
    for (testcase = 0; testcase < NOSUCH_TEST; testcase++) {
      /* Loop over radius */
      for (irad = 0; irad < nrad; irad++) {
        size_t rfail;
        /* Calculate radius */
        r = exp(h * irad);
        /* Compute the input data */
        compute_input(&values, &nspin, testcase, r);
        /* Compute the functional */
        rfail = check_xc(id, nspin, values, threshold);
        nfail += rfail;
#ifdef XC_DEBUG
        if (rfail) {
          printf("%lu failures for testcase = %i at r = %e\n", rfail,
                 (long unsigned)testcase, r);
        }
#endif
      }
    }

    if (nfail == 0) {
      /* Values are still OK, go further */
      threshold /= 10.0;
    } else {
      /* Working threshold */
      double ok_thresh = threshold*10.0;
      /* Current threshold */
      double cur_thresh = get_threshold(id);
      
#ifdef XC_DEBUG
      printf("%lu non-finite values encountered for threshold=%e\n",
             (long unsigned)nfail, threshold);
#endif

      printf("Current default threshold: %e\n", cur_thresh);
      /* We got numerical instabilities, so the threshold is the previous one */
      printf("Estimated working threshold: %e\n", ok_thresh);

      if(cur_thresh < ok_thresh) {
        printf("WARNING: FUNCTIONAL THRESHOLD FOR %s (ID = %i) SHOULD BE INCREASED TO %e!\n",fname,id,ok_thresh);
      }
      if(cur_thresh > ok_thresh) {
        printf("FUNCTIONAL THRESHOLD FOR %s (ID = %i) COULD BE DECREASED FROM %e TO %e!\n",fname,id,cur_thresh,ok_thresh);
      }
      
      break;
    }
  }

  libxc_free(fname);

  return 0;
}
