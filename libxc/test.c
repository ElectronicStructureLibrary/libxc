#include <math.h>
#include <stdio.h>

#include "util.h"
#include "xc.h"

void test_lda()
{
  lda_type l1, l2, l3;
  double rs;
  
  lda_init(&l1, XC_LDA_C_OB_PZ, XC_POLARIZED);
  lda_init(&l2, XC_LDA_C_GL, XC_POLARIZED);
  lda_init(&l3, XC_LDA_C_VWN, XC_POLARIZED);
  
  for(rs=0.001; rs<10; rs+=0.1){
    double ec, vc[2], rho[2];
    
    rho[0] = 1.0/(4.0/3.0*M_PI*pow(rs,3));
    
    rho[1] = 0.9*rho[0];
    rho[0] = 0.1*rho[0];
    
    lda(&l1, rho, &ec, vc);
    printf("%lf\t%lf\t%lf", rs, ec, vc[1]);
    
    lda(&l2, rho, &ec, vc);
    printf("\t%lf\t%lf", ec, vc[1]);
    
    lda(&l3, rho, &ec, vc);
    printf("\t%lf\t%lf\n", ec, vc[1]);
    
    
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
  test_tpss();

  return 0;
}
