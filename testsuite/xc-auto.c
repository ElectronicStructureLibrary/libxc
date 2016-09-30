/*
  Copyright (C) 2006-2007 M.A.L. Marques
  Copyright (C) 2014 Susi Lehtola

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <xc.h>
#include <xc_config.h>

/* Buffer size (line length) for file reads */
#define BUFSIZE 1024

typedef struct {
  /* Amount of data points */
  int n;

  /* Input: density, gradient, laplacian and kinetic energy density */
  FLOAT *rho;
  FLOAT *sigma;
  FLOAT *lapl;
  FLOAT *tau;

  /* Output: energy density and potentials for density, gradient, laplacian and tau */
  FLOAT *zk;
  FLOAT *vrho;
  FLOAT *vsigma;
  FLOAT *vlapl;
  FLOAT *vtau;

  /* ... and second derivatives */
  FLOAT *v2rho2;
  FLOAT *v2tau2;
  FLOAT *v2lapl2;
  FLOAT *v2rhotau;
  FLOAT *v2rholapl;
  FLOAT *v2lapltau;
  FLOAT *v2sigma2;
  FLOAT *v2rhosigma;
  FLOAT *v2sigmatau;
  FLOAT *v2sigmalapl;

  /* ... and third derivatives */
  FLOAT *v3rho3;
  FLOAT *v3rho2sigma;
  FLOAT *v3rhosigma2;
  FLOAT *v3sigma3;
} values_t;

void allocate_memory(values_t *data, int nspin) {
  data->rho         = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  data->sigma       = calloc(((nspin == 0) ? 1 : 3)*data->n, sizeof(FLOAT));
  data->lapl        = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  data->tau         = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  
  data->zk          = calloc(((nspin == 0) ? 1 : 1)*data->n, sizeof(FLOAT));
  data->vrho        = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  data->vsigma      = calloc(((nspin == 0) ? 1 : 3)*data->n, sizeof(FLOAT));
  data->vlapl       = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  data->vtau        = calloc(((nspin == 0) ? 1 : 2)*data->n, sizeof(FLOAT));
  
  data->v2rho2      = calloc(((nspin == 0) ? 1 : 3)*data->n, sizeof(FLOAT));
  data->v2tau2      = calloc(((nspin == 0) ? 1 : 3)*data->n, sizeof(FLOAT));
  data->v2lapl2     = calloc(((nspin == 0) ? 1 : 3)*data->n, sizeof(FLOAT));
  data->v2rhotau    = calloc(((nspin == 0) ? 1 : 4)*data->n, sizeof(FLOAT));
  data->v2rholapl   = calloc(((nspin == 0) ? 1 : 4)*data->n, sizeof(FLOAT));
  data->v2lapltau   = calloc(((nspin == 0) ? 1 : 6)*data->n, sizeof(FLOAT));
  data->v2sigma2    = calloc(((nspin == 0) ? 1 : 6)*data->n, sizeof(FLOAT));
  data->v2rhosigma  = calloc(((nspin == 0) ? 1 : 6)*data->n, sizeof(FLOAT));
  data->v2sigmatau  = calloc(((nspin == 0) ? 1 : 6)*data->n, sizeof(FLOAT));
  data->v2sigmalapl = calloc(((nspin == 0) ? 1 : 6)*data->n, sizeof(FLOAT));
  
  data->v3rho3      = calloc(((nspin == 0) ? 1 : 4)*data->n, sizeof(FLOAT));
  data->v3rho2sigma = calloc(((nspin == 0) ? 1 : 9)*data->n, sizeof(FLOAT));
  data->v3rhosigma2 = calloc(((nspin == 0) ? 1 :12)*data->n, sizeof(FLOAT));
  data->v3sigma3    = calloc(((nspin == 0) ? 1 :10)*data->n, sizeof(FLOAT));
}

void free_memory(values_t val) {
  free(val.rho); free(val.sigma); free(val.lapl); free(val.tau);

  free(val.zk); free(val.vrho); free(val.vsigma); free(val.vlapl); free(val.vtau);

  free(val.v2rho2); free(val.v2tau2); free(val.v2lapl2); free(val.v2rhotau);
  free(val.v2rholapl); free(val.v2lapltau); free(val.v2sigma2); free(val.v2rhosigma);
  free(val.v2sigmatau); free(val.v2sigmalapl);

  free(val.v3rho3); free(val.v3rho2sigma); free(val.v3rhosigma2); free(val.v3sigma3);
}

void read_data(const char *file, int nspin, values_t *data) {
  /* Format string */
#ifdef SINGLE_PRECISION
  static const char fmt[]="%f %f %f %f %f %f %f %f %f";
#else
  static const char fmt[]="%lf %lf %lf %lf %lf %lf %lf %lf %lf";
#endif

  /* Data buffer */
  char buf[BUFSIZE];
  char *cp;
  /* Input data file */
  FILE *in;
  /* Loop index */
  int i;
  /* Amount of points succesfully read */
  int nsucc;

  /* Helper variables */
  FLOAT rhoa, rhob;
  FLOAT sigmaaa, sigmaab, sigmabb;
  FLOAT lapla, laplb;
  FLOAT taua, taub;

  /* Open file */
  in=fopen(file,"r");
  if(!in) {
    fprintf(stderr,"Error opening input file %s.\n",file);
    exit(3);
  }

  /* Read amount of data points */
  cp = fgets(buf, BUFSIZE, in);
  if(cp != buf) {
    fprintf(stderr, "Error reading amount of data points.\n");
    exit(5);
  }
  nsucc = sscanf(buf, "%i", &(data->n));
  if(nsucc != 1) {
    fprintf(stderr, "Error reading amount of input data points.\n");
    exit(4);
  }

  /* Allocate memory */
  allocate_memory(data, nspin);

  for(i=0; i<data->n; i++) {
    /* Next line of input */
    cp = fgets(buf, BUFSIZE, in);
    if(cp != buf) {
      fprintf(stderr,"Read error on line %i.\n",i+1);
      exit(5);
    }
    /* Read data */
    nsucc = sscanf(buf, fmt, &rhoa, &rhob, &sigmaaa, &sigmaab, &sigmabb,
                   &lapla, &laplb, &taua, &taub);

    /* Error control */
    if(nsucc != 9) {
      fprintf(stderr, "Read error on line %i: only %i entries read.\n",
              i+1, nsucc);
      exit(5);
    }

    /* Store data (if clause suboptimal here but better for code clarity) */
    if(nspin==XC_POLARIZED) {
      data->rho[2*i]     = rhoa;
      data->rho[2*i+1]   = rhob;
      data->sigma[3*i]   = sigmaaa;
      data->sigma[3*i+1] = sigmaab;
      data->sigma[3*i+2] = sigmabb;
      data->lapl[2*i]    = lapla;
      data->lapl[2*i+1]  = laplb;
      data->tau[2*i]     = taua;
      data->tau[2*i+1]   = taub;
    } else {
      /* Construct full density data from alpha and beta channels */
      data->rho[i]   = rhoa + rhob;
      data->sigma[i] = sigmaaa + sigmabb + 2.0*sigmaab;
      data->lapl[i]  = lapla + laplb;
      data->tau[i]   = taua + taub;
    }
  }

  /* Close input file */
  fclose(in);
}

/*----------------------------------------------------------*/
#define MAX_DIFF 1e-8
void compare(char *what, int nspin, int n, int channels, double *v1, double *v2)
{
  int ii, max_ii;
  double max_diff, new_diff;

  max_diff = 0.0;
  max_ii   = 0;

  for(ii=0; ii<channels*n; ii++){
    if(fabs(v1[ii]) < MAX_DIFF)
      new_diff = fabs(v1[ii] - v2[ii]);
    else
      new_diff = fabs((v1[ii] - v2[ii])/v1[ii]);

    if(new_diff > max_diff){
      max_diff = new_diff;
      max_ii   = ii;
    }
  }

  if(max_diff > MAX_DIFF){
    printf("%s (nspin %d, channel %d)\n", what, nspin, ii % channels);
    printf("  point %d: %14.10e <> %14.10e\n", max_ii/channels, v1[max_ii], v2[max_ii]);
  }
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int func_id, nspin, i;
  
  if(argc != 3) {
    fprintf(stderr, "Usage:\n%s funct input\n", argv[0]);
    exit(1);
  }

  /* Get functional id */
  func_id = XC(functional_get_number)(argv[1]);
  if(func_id <= 0) {
    fprintf(stderr, "Functional '%s' not found\n", argv[1]);
    exit(1);
  }

  for(nspin=1; nspin<=2; nspin++) {
    values_t d[2];
    XC(func_type) func;
    int flags, family, ii, max_ii, channels;
    double max_diff;

    /* Helpers for properties that may not have been implemented */
    FLOAT *zk, *vrho, *v2rho2, *v3rho3;

    /* Read in data */
    read_data(argv[2], nspin, &(d[0]));
    read_data(argv[2], nspin, &(d[1]));

    /* Initialize functional */
    if(XC(func_init)(&func, func_id, nspin)) {
      fprintf(stderr, "Functional '%d' (%s) not found.\n"
              "Please report a bug against functional_get_number.\n", func_id, argv[1]);
      exit(1);
    }
    /* Get flags */
    flags  = func.info->flags;
    family = func.info->family;

    for(ii=0; ii<2; ii++){
      XC(func_set_derivatives) (&func, ii+1);
      
      /* Set helpers */
      zk     = (flags & XC_FLAGS_HAVE_EXC) ? d[ii].zk     : NULL;
      vrho   = (flags & XC_FLAGS_HAVE_VXC) ? d[ii].vrho   : NULL;
      v2rho2 = (flags & XC_FLAGS_HAVE_FXC) ? d[ii].v2rho2 : NULL;
      v3rho3 = (flags & XC_FLAGS_HAVE_KXC) ? d[ii].v3rho3 : NULL;

      /* Evaluate xc functional */
      switch(family) {
      case XC_FAMILY_LDA:
        XC(lda)(&func, d[ii].n, d[ii].rho, zk, vrho, v2rho2, v3rho3);
        break;
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA:
        XC(gga)(&func, d[ii].n, d[ii].rho, d[ii].sigma, zk, vrho, d[ii].vsigma,
                v2rho2, d[ii].v2rhosigma, d[ii].v2sigma2, 
                v3rho3, d[ii].v3rho2sigma, d[ii].v3rhosigma2, d[ii].v3sigma3);
        break;
      case XC_FAMILY_MGGA:
      case XC_FAMILY_HYB_MGGA:
        XC(mgga)(&func, d[ii].n, d[ii].rho, d[ii].sigma, d[ii].lapl, d[ii].tau, 
                 zk, vrho, d[ii].vsigma, d[ii].vlapl, d[ii].vtau,
                 v2rho2, d[ii].v2sigma2, d[ii].v2lapl2, d[ii].v2tau2, d[ii].v2rhosigma, 
                 d[ii].v2rholapl, d[ii].v2rhotau, d[ii].v2sigmalapl, d[ii].v2sigmatau, d[ii].v2lapltau);
        break;
        
      default:
        fprintf(stderr,"Support for family %i not implemented.\n",family);
        exit(1);
      }
    }

    /* energy */
    if(flags & XC_FLAGS_HAVE_EXC)
      compare("energy", nspin, d[0].n, 1, d[0].zk, d[1].zk);
 
    if(flags & XC_FLAGS_HAVE_VXC){
      compare("vrho", nspin, d[0].n, (nspin == 0) ? 1 : 2, d[0].vrho, d[1].vrho);

      if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA))
        compare("vsigma", nspin, d[0].n, (nspin == 0) ? 1 : 3, d[0].vsigma, d[1].vsigma);
    }

    if(flags & XC_FLAGS_HAVE_VXC){
      compare("v2rho2", nspin, d[0].n, (nspin == 0) ? 1 : 3, d[0].v2rho2, d[1].v2rho2);

      if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)){
        compare("v2sigma2",   nspin, d[0].n, (nspin == 0) ? 1 : 6, d[0].v2sigma2,   d[1].v2sigma2);
        compare("v2rhosigma", nspin, d[0].n, (nspin == 0) ? 1 : 6, d[0].v2rhosigma, d[1].v2rhosigma);
      }
    }

    if(flags & XC_FLAGS_HAVE_FXC){
      compare("v3rho3", nspin, d[0].n, (nspin == 0) ? 1 : 4, d[0].v3rho3, d[1].v3rho3);

      if(family & (XC_FAMILY_GGA | XC_FAMILY_HYB_GGA | XC_FAMILY_MGGA | XC_FAMILY_HYB_MGGA)){
        compare("v3rho2sigma", nspin, d[0].n, (nspin == 0) ? 1 : 9, d[0].v3rho2sigma, d[1].v3rho2sigma);
        compare("v3rhosigma2", nspin, d[0].n, (nspin == 0) ? 1 :12, d[0].v3rhosigma2, d[1].v3rhosigma2);
        compare("v3sigma3",    nspin, d[0].n, (nspin == 0) ? 1 :10, d[0].v3sigma3,    d[1].v3sigma3);
      }
    }
   
    XC(func_end)(&func);
    free_memory(d[0]);
    free_memory(d[1]);
  }
  
  return 0;
}
