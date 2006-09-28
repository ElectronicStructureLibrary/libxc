#include <math.h>
#include <stdio.h>

#include "../xc.h"

typedef struct {
  int functional;
  int nspin;

  double rho[2];        /* rhoa, rhob */
  double sigma[3];      /* sigmaaa, sigmaab, sigmabb */
  double zk;            /* energy density per unit particle */
  double vrho[2];       /* vrhoa, vrhob */
  double vsigma[3];     /* vsigmaaa, vsigmaab, vsigmabb */
  double v2rho[3];      /* v2rhoa2, v2rhoab, v2rhob2 */
  double v2rhosigma[6]; /* v2rhoasigmaaa, v2rhoasigmaab, v2rhoasigmabb
			   v2rhobsigmaaa, v2rhobsigmaab, v2rhobsigmabb */
  double v2sigma[6];    /* v2sigmaaa2, v2sigmaaaab, v2sigmaaabb
			   v2sigmaab2, v2sigmaabbb, v2sigmabb2 */
} xc_values_type;

/*----------------------------------------------------------*/
void init_values(xc_values_type *xc_values, char *argv[])
{
  xc_values->functional = atoi(argv[1]);
  xc_values->nspin      = atoi(argv[2]);
  xc_values->rho[0]     = atof(argv[3]);
  xc_values->rho[1]     = atof(argv[4]);
  xc_values->sigma[0]   = atof(argv[5]);
  xc_values->sigma[1]   = atof(argv[6]);
  xc_values->sigma[2]   = atof(argv[7]);

  int i;
  xc_values->zk = 0;
  for(i=0; i<2; i++){
    xc_values->vrho[i] = 0;
  }
  for(i=0; i<3; i++){
    xc_values->vsigma[i] = 0;
    xc_values->v2rho[i]  = 0;
  }

  for(i=0; i<6; i++){
    xc_values->v2rhosigma[i]  = 0;
    xc_values->v2sigma[i]     = 0;
  }
}


/*----------------------------------------------------------*/
void print_values(xc_values_type *xc)
{
  printf(" rhoa= %#0.2E rhob= %#0.2E sigmaaa= %#0.2E sigmaab= %#0.2E sigmabb= %#0.2E\n",
	 xc->rho[0], xc->rho[1],
	 xc->sigma[0], xc->sigma[1], xc->sigma[2]);
  printf(" zk            = %#19.12E\n"
	 " vrhoa         = %#19.12E\n"
	 " vrhob         = %#19.12E\n"
	 " vsigmaaa      = %#19.12E\n"
	 " vsigmaab      = %#19.12E\n"
	 " vsigmabb      = %#19.12E\n\n"
	 " v2rhoa2       = %#19.12E\n"
	 " v2rhoab       = %#19.12E\n"
	 " v2rhob2       = %#19.12E\n\n"
	 " v2rhoasigmaaa = %#19.12E\n"
	 " v2rhoasigmaab = %#19.12E\n"
	 " v2rhoasigmabb = %#19.12E\n"
	 " v2rhobsigmaaa = %#19.12E\n"
	 " v2rhobsigmaab = %#19.12E\n"
	 " v2rhobsigmabb = %#19.12E\n\n"
	 " v2sigmaaa2    = %#19.12E\n"
	 " v2sigmaaaab   = %#19.12E\n"
	 " v2sigmaaabb   = %#19.12E\n"
	 " v2sigmaab2    = %#19.12E\n"
	 " v2sigmaabbb   = %#19.12E\n"
	 " v2sigmabb2    = %#19.12E\n",
	 xc->zk, 
	 xc->vrho[0], xc->vrho[1],
	 xc->vsigma[0], xc->vsigma[1], xc->vsigma[2],
	 xc->v2rho[0], xc->v2rho[1], xc->v2rho[2],
	 xc->v2rhosigma[0], xc->v2rhosigma[1], xc->v2rhosigma[2],
	 xc->v2rhosigma[3], xc->v2rhosigma[4], xc->v2rhosigma[5],
	 xc->v2sigma[0], xc->v2sigma[1], xc->v2sigma[2],
	 xc->v2sigma[3], xc->v2sigma[4], xc->v2sigma[5]
	 );
}


/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  xc_values_type xc;
  lda_type lda_func;
  gga_type gga_func;

  if(argc != 8){
    printf("Usage:\n%s funct pol rhoa rhob sigmaaa sigmaab sigmabb\n", argv[0]);
    return 1;
  }

  init_values(&xc, argv);

  if(xc.nspin == 1){
    xc.rho[0]    += xc.rho[1];
    xc.sigma[0] += 2.0*xc.sigma[1] + xc.sigma[2];
  }

  if(xc.functional < 100){
    if(xc.functional == XC_LDA_X){
      lda_x_init(&lda_func, xc.nspin, 3, 0);
    }else{
      lda_init(&lda_func, xc.functional, xc.nspin);
    }

    lda(&lda_func, xc.rho, &xc.zk, xc.vrho, NULL);

  }else if(xc.functional < 200){
    gga_init(&gga_func, xc.functional, xc.nspin);
    gga(&gga_func, xc.rho, xc.sigma, &xc.zk, xc.vrho, xc.vsigma);
    gga_end(&gga_func);
  }

  if(xc.nspin == 1){
    xc.vrho[1] = xc.vrho[0];
    xc.zk *= xc.rho[0] ;
  }else{
    xc.zk *= (xc.rho[0] + xc.rho[1]);
  }

  print_values(&xc);

  return 0;
}

