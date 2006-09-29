#include <math.h>
#include <stdio.h>

#include "util.h"
#include "xc.h"

void test_lda()
{
  lda_type l1, l2, l3;
  double rs;
  
  lda_x_init(&l1, XC_POLARIZED, 3);
  lda_init(&l2, XC_LDA_C_PZ, XC_POLARIZED);
  lda_init(&l3, XC_LDA_C_VWN, XC_POLARIZED);

  for(rs=10; rs>=0.1; rs-=0.01){
    double dens, zeta, rho[2];
    double ec1, vc1[2], fxc1[4];
    double ec2, vc2[2], fxc2[4];
    
    dens = 1.0/(4.0/3.0*M_PI*pow(rs,3));
    //zeta = 1;
    //rho[0] = dens*(1.0 + zeta)/2.0;
    //rho[1] = dens*(1.0 - zeta)/2.0;

    rho[0] = dens;
    rho[1] = 0;
    dens   = (rho[0] + rho[1]);
    zeta   = (rho[0] - rho[1])/dens;

    lda(&l2, rho, &ec1, vc1);
    lda(&l3, rho, &ec2, vc2);
    lda_fxc(&l1, rho, fxc1);
    lda_fxc(&l3, rho, fxc2);

    printf("%e\t%e\t%e\t%e\t%e\n", dens, fxc1[0], fxc1[1], fxc1[2], fxc1[3]);
    
  }
}

void test_tpss()
{
  mgga_type tpss;
  int i;

  mgga_init(&tpss, XC_MGGA_C_TPSS, XC_UNPOLARIZED);
  
  for(i=0; i<100; i++){
    double n, gr[3], tau;
    double e, dedd, dedgd[3], dedtau;

    n = 0.1;
    gr[0] = 0.01;
    gr[1] = 0.1;
    gr[2] = 0.1;
    tau   = 0.01 + i/100.0;

    mgga(&tpss, &n, gr, &tau, &e, &dedd, dedgd, &dedtau);
    printf("%16.10lf\t%16.10lf\t%16.10lf\n", tau, n*e, dedtau);
  }
}

int main()
{
  test_lda();

  return 0;
}
