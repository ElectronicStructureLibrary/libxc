/*
 Copyright (C) 2006-2007 M.A.L. Marques
 Copyright (C) 2014 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <xc.h>

/* Buffer size (line length) for file reads */
#define BUFSIZE 1024

typedef struct {
  /* Amount of data points */
  int n;

  /* Input: density, gradient, laplacian and kinetic energy density */
  double *rho;
  double *sigma;
  double *lapl;
  double *tau;

  /* Output: energy density */
  double *zk;

  /* .. and potentials for density, gradient, laplacian and tau */
  double *vrho;
  double *vsigma;
  double *vlapl;
  double *vtau;

  /* ... and second derivatives */
  double *v2rho2;
  double *v2tau2;
  double *v2lapl2;
  double *v2rhotau;
  double *v2rholapl;
  double *v2lapltau;
  double *v2sigma2;
  double *v2rhosigma;
  double *v2sigmatau;
  double *v2sigmalapl;

  /* ... and third derivatives */
  double *v3rho3;
} values_t;

void allocate_memory(values_t *data, int nspin, int order)
{
  data->zk = NULL;
  data->vrho = NULL;
  data->vsigma = NULL;
  data->vlapl = NULL;
  data->vtau = NULL;
  data->v2rho2 = NULL;
  data->v2tau2 = NULL;
  data->v2lapl2 = NULL;
  data->v2rhotau = NULL;
  data->v2rholapl = NULL;
  data->v2lapltau = NULL;
  data->v2sigma2 = NULL;
  data->v2rhosigma = NULL;
  data->v2sigmatau = NULL;
  data->v2sigmalapl = NULL;
  data->v3rho3 = NULL;

  switch(nspin) {
    case (XC_UNPOLARIZED):
      data->rho = calloc(data->n, sizeof(double));
      data->sigma = calloc(data->n, sizeof(double));
      data->lapl = calloc(data->n, sizeof(double));
      data->tau = calloc(data->n, sizeof(double));
      switch (order) {
        case (0):
          data->zk = calloc(data->n, sizeof(double));
          break;
        case (1):
          data->vrho = calloc(data->n, sizeof(double));
          data->vsigma = calloc(data->n, sizeof(double));
          data->vlapl = calloc(data->n, sizeof(double));
          data->vtau = calloc(data->n, sizeof(double));
          break;
        case (2):
          data->v2rho2 = calloc(data->n, sizeof(double));
          data->v2tau2 = calloc(data->n, sizeof(double));
          data->v2lapl2 = calloc(data->n, sizeof(double));
          data->v2rhotau = calloc(data->n, sizeof(double));
          data->v2rholapl = calloc(data->n, sizeof(double));
          data->v2lapltau = calloc(data->n, sizeof(double));
          data->v2sigma2 = calloc(data->n, sizeof(double));
          data->v2rhosigma = calloc(data->n, sizeof(double));
          data->v2sigmatau = calloc(data->n, sizeof(double));
          data->v2sigmalapl = calloc(data->n, sizeof(double));
          break;
        case (3):
          data->v3rho3 = calloc(data->n, sizeof(double));
          break;
        default:
          fprintf(stderr, "order = %i not recognized.\n", order);
          exit(2);
      }
      break;

    case (XC_POLARIZED):
      data->rho = calloc(2*data->n, sizeof(double));
      data->sigma = calloc(3*data->n, sizeof(double));
      data->lapl = calloc(2*data->n, sizeof(double));
      data->tau = calloc(2*data->n, sizeof(double));
      switch (order) {
        case (0):
          data->zk = calloc(data->n, sizeof(double));
          break;
        case (1):
          data->vrho = calloc(2*data->n, sizeof(double));
          data->vsigma = calloc(3*data->n, sizeof(double));
          data->vlapl = calloc(2*data->n, sizeof(double));
          data->vtau = calloc(2*data->n, sizeof(double));
          break;
        case (2):
          data->v2rho2 = calloc(3*data->n, sizeof(double));
          data->v2tau2 = calloc(3*data->n, sizeof(double));
          data->v2lapl2 = calloc(3*data->n, sizeof(double));
          data->v2rhotau = calloc(4*data->n, sizeof(double));
          data->v2rholapl = calloc(4*data->n, sizeof(double));
          data->v2lapltau = calloc(4*data->n, sizeof(double));
          data->v2sigma2 = calloc(6*data->n, sizeof(double));
          data->v2rhosigma = calloc(6*data->n, sizeof(double));
          data->v2sigmatau = calloc(6*data->n, sizeof(double));
          data->v2sigmalapl = calloc(6*data->n, sizeof(double));
          break;
        case (3):
          data->v3rho3 = calloc(4*data->n, sizeof(double));
          break;
        default:
          fprintf(stderr, "order = %i not recognized.\n", order);
          exit(2);
      }
      break;

    default:
      fprintf(stderr, "nspin = %i not recognized.\n", nspin);
      exit(2);
  }

}

void free_memory(values_t val)
{
  free(val.rho);
  free(val.sigma);
  free(val.lapl);
  free(val.tau);
  free(val.zk);
  free(val.vrho);
  free(val.vsigma);
  free(val.vlapl);
  free(val.vtau);
  free(val.v2rho2);
  free(val.v2tau2);
  free(val.v2lapl2);
  free(val.v2rhotau);
  free(val.v2rholapl);
  free(val.v2lapltau);
  free(val.v2sigma2);
  free(val.v2rhosigma);
  free(val.v2sigmatau);
  free(val.v2sigmalapl);
  free(val.v3rho3);
}

values_t read_data(const char *file, int nspin, int order) {
  /* Format string */
  static const char fmt[]="%lf %lf %lf %lf %lf %lf %lf %lf %lf";

  /* Data buffer */
  char buf[BUFSIZE];
  char *cp;
  /* Input data file */
  FILE *in;
  /* Loop index */
  int i;
  /* Amount of points succesfully read */
  int nsucc;
  /* Returned data */
  values_t data;

  /* Helper variables */
  double rhoa, rhob;
  double sigmaaa, sigmaab, sigmabb;
  double lapla, laplb;
  double taua, taub;

  /* Open file */
  in=fopen(file,"r");
  if(!in) {
    fprintf(stderr,"Error opening input file %s.\n",file);
    exit(3);
  }

  /* Read amount of data points */
  cp=fgets(buf,BUFSIZE,in);
  if(cp!=buf) {
    fprintf(stderr,"Error reading amount of data points.\n");
    exit(5);
  }
  nsucc=sscanf(buf,"%i",&data.n);
  if(nsucc!=1) {
    fprintf(stderr,"Error reading amount of input data points.\n");
    exit(4);
  }

  /* Allocate memory */
  allocate_memory(&data, nspin, order);

  for(i=0;i<data.n;i++) {
    /* Next line of input */
    cp=fgets(buf,BUFSIZE,in);
    if(cp!=buf) {
      fprintf(stderr,"Read error on line %i.\n",i+1);
      free_memory(data);
      exit(5);
    }
    /* Read data */
    nsucc=sscanf(buf, fmt, &rhoa, &rhob, &sigmaaa, &sigmaab, &sigmabb,	\
		 &lapla, &laplb, &taua, &taub);

    /* Error control */
    if(nsucc!=9) {
      fprintf(stderr,"Read error on line %i: only %i entries read.\n",i+1,nsucc);
      free_memory(data);
      exit(5);
    }

    /* Store data (if clause suboptimal here but better for code clarity) */
    if(nspin==XC_POLARIZED) {
      data.rho[2*i]=rhoa;
      data.rho[2*i+1]=rhob;
      data.sigma[3*i]=sigmaaa;
      data.sigma[3*i+1]=sigmaab;
      data.sigma[3*i+2]=sigmabb;
      data.lapl[2*i]=lapla;
      data.lapl[2*i+1]=laplb;
      data.tau[2*i]=taua;
      data.tau[2*i+1]=taub;
    } else {
      /* Construct full density data from alpha and beta channels */
      data.rho[i]=rhoa + rhob;
      data.sigma[i]=sigmaaa + sigmabb + 2.0*sigmaab;
      data.lapl[i]=lapla + laplb;
      data.tau[i]=taua + taub;
    }
  }

  /* Close input file */
  fclose(in);

  return data;
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int func_id, nspin, order, i;
  /* Helpers for properties that may not have been implemented */
  double *zk, *vrho, *v2rho2, *v3rho3;

  static const char efmt[] =" % .16e";
  static const char efmt2[]=" % .16e % .16e";
  static const char efmt3[]=" % .16e % .16e % .16e";
  static const char sfmt[] =" %23s";
  static const char sfmt2[]=" %23s %23s";
  static const char sfmt3[]=" %23s %23s %23s";

  if(argc != 6) {
    fprintf(stderr, "Usage:\n%s funct nspin order input output\n", argv[0]);
    exit(1);
  }

  /* Get functional id */
  func_id = xc_functional_get_number(argv[1]);
  if(func_id <= 0) {
    fprintf(stderr, "Functional '%s' not found\n", argv[1]);
    exit(1);
  }

  /* Spin-polarized or unpolarized ? */
  nspin = atoi(argv[2]);

  /* Order of derivatives to compute */
  order = atoi(argv[3]);

  /* Data array */
    values_t d;
  /* Functional evaluator */
  xc_func_type func;
  /* Flags for functional */
  int flags;
  /* Functional family */
  int family;
  /* Output file */
  FILE *out;
  /* Output file name */
  char *fname;

  /* Read in data */
  d = read_data(argv[4], nspin, order);

  /* Initialize functional */
  if(xc_func_init(&func, func_id, nspin)) {
    fprintf(stderr, "Functional '%d' (%s) not found.\nPlease report a bug against functional_get_number.\n", func_id, argv[1]);
    exit(1);
  }
  /* Get flags */
  flags=func.info->flags;
  family=func.info->family;

  /* Set helpers */
  zk     = (flags & XC_FLAGS_HAVE_EXC) ? d.zk     : NULL;
  vrho   = (flags & XC_FLAGS_HAVE_VXC) ? d.vrho   : NULL;
  v2rho2 = (flags & XC_FLAGS_HAVE_FXC) ? d.v2rho2 : NULL;
  v3rho3 = (flags & XC_FLAGS_HAVE_KXC) ? d.v3rho3 : NULL;

  /* Evaluate xc functional */
  switch(family) {
  case XC_FAMILY_LDA:
    xc_lda(&func, d.n, d.rho, zk, vrho, v2rho2, v3rho3);
    break;
  case XC_FAMILY_GGA:
  case XC_FAMILY_HYB_GGA:
    xc_gga(&func, d.n, d.rho, d.sigma, zk, vrho, d.vsigma,		\
    v2rho2, d.v2rhosigma, d.v2sigma2, NULL, NULL, NULL, NULL);
    break;
  case XC_FAMILY_MGGA:
  case XC_FAMILY_HYB_MGGA:
    xc_mgga(&func, d.n, d.rho, d.sigma, d.lapl, d.tau, zk, vrho, d.vsigma, d.vlapl, d.vtau, \
     v2rho2, d.v2sigma2, d.v2lapl2, d.v2tau2, d.v2rhosigma, d.v2rholapl, d.v2rhotau, \
     d.v2sigmalapl, d.v2sigmatau, d.v2lapltau);
    break;

  default:
    fprintf(stderr,"Support for family %i not implemented.\n",family);
    free_memory(d);
    exit(1);
  }

  /* Open output file */
  fname = argv[5];
  out = fopen(fname,"w");
  if(!out) {
    fprintf(stderr,"Error opening output file %s.\n",fname);
    free_memory(d);
    exit(1);
  }

  /* Functional id and amount of lines in output */
  fprintf(out, "%i %i %i\n", func_id, d.n, order);

  switch (order) {
    case (0): /* energy */
      fprintf(out, sfmt, "zk");
      break;
    case (1): /* first order derivatives */
      if (nspin == XC_POLARIZED) {
        fprintf(out, sfmt2, "vrho(a)", "vrho(b)");
        if (family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA))
          fprintf(out, sfmt3, "vsigma(aa)", "vsigma(ab)", "vsigma(bb)");
        if (family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt2, "vlapl(a)", "vlapl(b)");
          fprintf(out, sfmt2, "vtau(a)", "vtau(b)");
        }
      } else {
        fprintf(out, sfmt, "vrho");
        if (family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA))
          fprintf(out, sfmt, "vsigma");
        if(family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt, "vlapl");
          fprintf(out, sfmt, "vtau");
        }
      }
      break;

    case (2): /* second order derivatives */
      if (nspin == XC_POLARIZED) {
        fprintf(out,sfmt3,"v2rho(aa)","v2rho(ab)","v2rho(bb)");
        if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt3, "v2sigma2(aa-aa)", "v2sigma2(aa-ab)", "v2sigma2(aa-bb)");
          fprintf(out, sfmt3, "v2sigma2(ab-ab)", "v2sigma2(ab-bb)", "v2sigma2(bb-bb)");
          fprintf(out, sfmt3, "v2rho(a)sigma(aa)", "v2rho(a)sigma(ab)", "v2rho(a)sigma(bb)");
          fprintf(out, sfmt3, "v2rho(b)sigma(aa)", "v2rho(b)sigma(ab)", "v2rho(b)sigma(bb)");
        }
        if(family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt3, "v2lapl2(aa)", "v2lapl2(ab)", "v2lapl2(bb)");
          fprintf(out, sfmt3, "v2tau2(aa)", "v2tau2(ab)", "v2tau2(bb)");
          fprintf(out, sfmt3, "v2rholapl(aa)", "v2rholapl(ab)", "v2rholapl(bb)");
          fprintf(out, sfmt3, "v2rhotau(aa)", "v2rhotau(ab)", "v2rhotau(bb)");
          fprintf(out, sfmt3, "v2lapltau(aa)", "v2lapltau(ab)", "v2lapltau(bb)");
          fprintf(out, sfmt3, "v2sigma(aa)tau(a)", "v2sigma(aa)tau(b)", "v2sigma(ab)tau(a)");
          fprintf(out, sfmt3, "v2sigma(ab)tau(b)", "v2sigma(bb)tau(a)", "v2sigma(bb)tau(b)");
          fprintf(out, sfmt3, "v2sigma(aa)lapl(a)", "v2sigma(aa)lapl(b)", "v2sigma(ab)lapl(a)");
          fprintf(out, sfmt3, "v2sigma(ab)lapl(b)", "v2sigma(bb)lapl(a)", "v2sigma(bb)lapl(b)");
        }
      } else {
        fprintf(out,sfmt,"v2rho");
        if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt, "v2sigma2");
          fprintf(out, sfmt, "v2rhosigma");
        }

        if(family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
          fprintf(out, sfmt, "v2lapl2");
          fprintf(out, sfmt, "v2tau2");
          fprintf(out, sfmt, "v2rholapl");
          fprintf(out, sfmt, "v2rhotau");
          fprintf(out, sfmt, "v2lapltau");
          fprintf(out, sfmt, "v2sigmatau");
          fprintf(out, sfmt, "v2sigmalapl");
        }
      }
      break;

    default: /* higher order derivatives ... to be done */
      fprintf(stderr, "order = %i not recognized.\n", order);
      exit(2);
  }
  fprintf(out,"\n");

  /* Loop over data points */
  for(i=0;i<d.n;i++) {

    switch (order) {
      case (0): /* energy */
        fprintf(out, efmt, d.zk[i]);
        break;
      case (1): /* first order derivatives */
        if (nspin == XC_POLARIZED) {
          fprintf(out, efmt2, d.vrho[2 * i], d.vrho[2 * i + 1]);
          if (family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA))
            fprintf(out, efmt3, d.vsigma[3 * i], d.vsigma[3 * i + 1], d.vsigma[3 * i + 2]);
          if (family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt2, d.vlapl[2 * i], d.vlapl[2 * i + 1]);
            fprintf(out, efmt2, d.vtau[2 * i], d.vtau[2 * i + 1]);
          }
        } else {
          fprintf(out, efmt, d.vrho[i]);
          if (family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA))
            fprintf(out, efmt, d.vsigma[i]);
          if (family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt, d.vlapl[i]);
            fprintf(out, efmt, d.vtau[i]);
          }
        }
        break;

      case (2): /* second order derivatives */
        if (nspin == XC_POLARIZED) {
          fprintf(out, efmt3, d.v2rho2[3*i], d.v2rho2[3*i + 1], d.v2rho2[3*i + 2]);
          if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt3, d.v2sigma2[6*i], d.v2sigma2[6*i + 1], d.v2sigma2[6*i + 2]);
            fprintf(out, efmt3, d.v2sigma2[6*i + 3], d.v2sigma2[6*i + 4], d.v2sigma2[6*i + 5]);
            fprintf(out, efmt3, d.v2rhosigma[6*i], d.v2rhosigma[6*i + 1], d.v2rhosigma[6*i + 2]);
            fprintf(out, efmt3, d.v2rhosigma[6*i + 3], d.v2rhosigma[6*i + 4], d.v2rhosigma[6*i + 5]);
          }
          if(family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt3, d.v2lapl2[3*i], d.v2lapl2[3*i + 1], d.v2lapl2[3*i + 2]);
            fprintf(out, efmt3, d.v2tau2[3*i], d.v2tau2[3*i + 1], d.v2tau2[3*i + 2]);
            fprintf(out, efmt3, d.v2rholapl[3*i], d.v2rholapl[3*i + 1], d.v2rholapl[3*i + 2]);
            fprintf(out, efmt3, d.v2rhotau[3*i], d.v2rhotau[3*i + 1], d.v2rhotau[3*i + 2]);
            fprintf(out, efmt3, d.v2lapltau[3*i], d.v2lapltau[3*i + 1], d.v2lapltau[3*i + 2]);
            fprintf(out, efmt3, d.v2sigmatau[3*i], d.v2sigmatau[3*i + 1], d.v2sigmatau[3*i + 2]);
            fprintf(out, efmt3, d.v2sigmatau[3*i + 3], d.v2sigmatau[3*i + 4], d.v2sigmatau[3*i + 5]);
            fprintf(out, efmt3, d.v2sigmalapl[3*i], d.v2sigmalapl[3*i + 1], d.v2sigmalapl[3*i + 2]);
            fprintf(out, efmt3, d.v2sigmalapl[3*i + 3], d.v2sigmalapl[3*i + 4], d.v2sigmalapl[3*i + 5]);
          }
        } else {
          fprintf(out, efmt, d.v2rho2[i]);
          if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt, d.v2sigma2[i]);
            fprintf(out, efmt, d.v2rhosigma[i]);
          }
          if(family & (XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)) {
            fprintf(out, efmt, d.v2lapl2[i]);
            fprintf(out, efmt, d.v2tau2[i]);
            fprintf(out, efmt, d.v2rholapl[i]);
            fprintf(out, efmt, d.v2rhotau[i]);
            fprintf(out, efmt, d.v2lapltau[i]);
            fprintf(out, efmt, d.v2sigmatau[i]);
            fprintf(out, efmt, d.v2sigmalapl[i]);
          }
        }
        break;

     default: /* higher order derivatives ... to be done */
        fprintf(stderr, "order = %i not recognized.\n", order);
        exit(2);
    }

    fprintf(out,"\n");
  }

  xc_func_end(&func);
  free_memory(d);
  fclose(out);

  return 0;
}
