#include <math.h>
#include <stdio.h>

#include "util.h"
#include "xc.h"

void test_lda()
{
  lda_type l1, l2, l3;
  double rs;
  
  lda_x_init(&l1, XC_POLARIZED, 3);
  lda_init(&l2, XC_LDA_C_PW, XC_POLARIZED);
  lda_init(&l3, XC_LDA_C_VWN, XC_POLARIZED);
  
  for(rs=0.01; rs<2; rs+=0.001){
    double dens, zeta, ec, vc[2], rho[2], fxc[4];
    
    //dens = 1.0/(4.0/3.0*M_PI*pow(rs,3));
    //zeta = 0.5;
    //rho[0] = dens*(1.0 + zeta)/2.0;
    //rho[1] = dens*(1.0 - zeta)/2.0;

    rho[0] = rs;
    rho[1] = 0.1;
    dens   = (rho[0] + rho[1]);
    zeta   = (rho[0] - rho[1])/dens;

    lda(&l2, rho, &ec, vc);
    lda_fxc(&l2, rho, fxc);

    printf("%e\t%e\t%e\t%e\n", rho[0], dens*ec, vc[0], fxc __(0,0));
    
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
