/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

#include "util.h"
#include "xc.h"

void func(double *x, int n, void *ex)
{
  int i;
  for(i=0; i<n;i++)
    x[i] = cos(x[i]);
}

void test_integration()
{
  double a, b, result;

  for(b=1e-8; b<5; b+=0.001){
    result = xc_integrate(func, NULL, a, b);
    printf("%lf %lf\n", b, result);
  }
}

#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_bessel.h>

void test_expi()
{
  double b, r1, r2, r3, r4;
  int n;

  n = 1;
  for(b=1e-3; b<50; b+=0.01){
    r1 = xc_bessel_K0(b);
    //r2 = gsl_sf_bessel_K0(b);
    //r3 = bessi1(b);
    printf("%5.3lf %12.10lf %12.10lf %12.10lf\n", b, r1, r2, r3);
  }  
}

void test_lda()
{
  XC(func_type) l1, l2, l3;
  int i;
  
  //XC(func_init)(&l1, XC_LDA_XC_KSDT2, XC_POLARIZED);
  //XC(func_init)(&l1, XC_LDA_C_PW, XC_UNPOLARIZED);
  XC(func_init)(&l3, XC_LDA_X, XC_UNPOLARIZED);

  //XC(lda_xc_ksdt_set_params)(&l1, 1000.0);
  //XC(lda_c_1d_csc_set_params)(&l2, 1, 1.0);

  for(i=0; i<=1000; i++){
    double dens, rs, zeta, rho[2];
    double ec1, vc1[2], fxc1[3], kxc1[4];
    double ec2, vc2[2], fxc2[3], kxc2[4];
    double ec3, vc3[2], fxc3[3], kxc3[4];
    
    //rs   = 0.5 + i/500.0;
    //zeta = -1.0 + 2.0*i/1000000000.0;

    //dens = 1.0/(4.0/3.0*M_PI*POW(rs,3)); /* 3D */
    //dens = 1.0/(2.0*rs); /* 1D */

    //rho[0] = dens*(1.0 + zeta)/2.0;
    //rho[1] = dens*(1.0 - zeta)/2.0;

    rho[0] = 0.05 + i/1000.0;
    rho[1] = 0.001;

    //rho[0] = 1.0/(2.0*rs);
    //rho[1] = 0.0;

    //dens = rho[0] + rho[1];

    XC(lda)(&l1, 1, rho, &ec1, vc1, NULL, NULL);
    //XC(lda)(&l2, 1, rho, &ec2, vc2, NULL, NULL);
    //XC(lda_fxc_fd)(&l2, rho, fxc2);
    //XC(lda_kxc_fd)(&l2, rho, kxc2);

    //rho[0] = dens; rho[1] = 0.0;
    //XC(lda)(&l3, rho, &ec3, vc3, fxc3, kxc3);

    // printf("%e\t%e\t%e\n", dens, (fxc1[0]+2.0*fxc1[1]+fxc1[2])/4.0, fxc3[0]);
    // printf("%e\t%e\t%e\n", dens, (kxc1[0]+3.0*kxc1[1]+3.0*kxc1[2]+kxc1[3])/8.0, kxc3[0]);

    printf("%e\t%e\t%e\n", rho[0], (rho[0] + rho[1])*ec1, vc1[0]);
  }
}

void test_ak13()
{
  XC(func_type) gga;
  double beta = 0.13;
  double x, rho[2] = {0.0, 0.0}, sigma[3] = {0.0, 0.0, 0.0}, zk, vrho[2], vsigma[3];
  double tmp1, tmp2;
  int i;

  XC(func_init)(&gga,  XC_GGA_X_AK13,  XC_POLARIZED);

  for(i=0; i<=10000; i++){

    x = 500.0*i/(10000.0);
    rho[0]   = 0.12*exp(-beta * x);
    sigma[0] = 0.12*0.12*beta*beta * rho[0]*rho[0];
    
    XC(gga_exc_vxc)(&gga,  1, rho, sigma, &zk, vrho, vsigma);

    tmp2 = 1.74959015598863046792081721182*beta*x/3.0- 1.62613336586517367779736042170*log(x);
    fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\t%16.10lf\n", x, vrho[0], vsigma[0]*sqrt(sigma[0]), 
	    -X_FACTOR_C*X2S*tmp2/2.0);
  }
  
}

void test_enhance()
{
  double x, f, dfdx, d2fdx2, d3fdx3;

  XC(func_type) gga1, gga2, gga3, gga4;

  XC(func_init)(&gga1,  XC_GGA_X_B88, XC_POLARIZED);
  XC(func_init)(&gga2, XC_GGA_X_GG99, XC_POLARIZED);

  for(x=0.01; x<200; x+=0.01){
    //printf("%le", x);
    //XC(gga_x_b88_enhance)(&gga1, 1, x, &f, &dfdx, &d2fdx2, &d3fdx3);
    //printf("\t%le", -f*X_FACTOR_C);
    XC(gga_x_gg99_enhance)(&gga2, 1, x, &f, &dfdx, &d2fdx2, &d3fdx3);
    printf("%le\t%le\t%le\n", x, f, dfdx);
  }
}


void test_gga()
{
  XC(func_type) gga, ggap;
  int i, npoints = 1;
  double *rho, *sigma;
  double *zk, zkp, *vrho, vrhop[2], *vsigma, vsigmap[3];
  double *v2rho2, *v2rhosigma, *v2sigma2;
  double *v3rho3, *v3rho2sigma, *v3rhosigma2, *v3sigma3;

  rho         = malloc( 2*npoints*sizeof(double));
  sigma       = malloc( 3*npoints*sizeof(double));
  zk          = malloc( 1*npoints*sizeof(double));
  vrho        = malloc( 2*npoints*sizeof(double));
  vsigma      = malloc( 3*npoints*sizeof(double));
  v2rho2      = malloc( 3*npoints*sizeof(double));
  v2rhosigma  = malloc( 6*npoints*sizeof(double));
  v2sigma2    = malloc( 6*npoints*sizeof(double));
  v3rho3      = malloc( 4*npoints*sizeof(double));
  v3rho2sigma = malloc( 9*npoints*sizeof(double));
  v3rhosigma2 = malloc(12*npoints*sizeof(double));
  v3sigma3    = malloc(10*npoints*sizeof(double));

  
  XC(func_init)(&gga,  XC_GGA_C_PW91,  XC_POLARIZED);

  /*
  for(i=1; i<=10000; i++){
    double x = 30.0*i/(10000.0), f, df, d2f, d3f;
    double c2 = 10.0/81.0;

    XC(gga_x_b88_enhance)(&gga, 3, x, &f, &df, &d2f, &d3f);
    printf("%20.14e %20.14e", X2S*X2S*x*x, (f-1.0)/(c2*X2S*X2S*x*x));
    XC(gga_x_ev93_enhance)(&gga, 3, x, &f, &df, &d2f, &d3f);
    printf(" %20.14e\n", (f-1.0)/(c2*X2S*X2S*x*x));
  }
  exit(0);
  */

   for(i=0; i<=10000; i++){
     rho[0]   = .01;
     rho[1]   = 0.2;
     sigma[0] = 0.1;
     sigma[1] = 0.00002;
     sigma[2] = 0.5 + i/1000.0;

     XC(gga)(&gga, 1, rho, sigma, zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, NULL, v3rho2sigma, v3rhosigma2, v3sigma3);

     fprintf(stderr, "%16.10le\t%16.10le\t%16.10le\n", sigma[2], vsigma[2], v2sigma2[5]);
   }

  /*
  for(i=0; i<npoints; i++){
    rho[2*i + 0]   = 0.048 + i/10000.0;
    rho[2*i + 1]   = 0.025;
    sigma[3*i + 0] = 0.0046;
    sigma[3*i + 1] = 0.0044;
    sigma[3*i + 2] = 0.0041;
  }

  XC(gga)(&gga,  npoints, rho, sigma, zk,  vrho,  vsigma,  v2rho2,  v2rhosigma,  v2sigma2, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3);

  for(i=0; i<npoints; i++){    
    fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", rho[2*i + 0], vrho[2*i + 0], v2rho2[3*i + 0]);
  }
  */

  XC(func_end)(&gga);

  free(rho); free(sigma);
  free(zk); free(vrho); free(vsigma);
  free(v2rho2); free(v2rhosigma); free(v2sigma2);
  free(v3rho3); free(v3rho2sigma); free(v3rhosigma2); free(v3sigma3);
}


void test_mgga()
{
  XC(func_type) mgga1, mgga2;
  int i;

  XC(func_init)(&mgga1, XC_GGA_C_LYP, XC_POLARIZED);
  XC(func_init)(&mgga2, XC_MGGA_X_SCAN, XC_POLARIZED);
  //XC(mgga_c_tpss_init)(tpss2.mgga);
  
  for(i=0; i<=1000; i++){
    double rho[2], sigma[3], tau[2], lapl[2];
    double zk,   vrho[2],  vsigma[3],  vtau[2],  vlapl[2];
    double zk2, vrho2[2], vsigma2[3], vtau2[2], vlapl2[2];
    double v2rho2[3], v2sigma2[6], v2lapl2[3], v2tau2[3];
    double v2rhosigma[6], v2rholapl[3], v2rhotau[3];
    double v2sigmalapl[6], v2sigmatau[6], v2lapltau[3];

    rho[0]   = 0.3;
    rho[1]   = 0.4;
    sigma[0] = 0.2032882206468622;
    sigma[1] = 0.11;
    sigma[2] = 0.7 + i/10.0;
    tau[0]   = 1.0;
    tau[1]   = 0.15;
    lapl[0]  = -0.1518421131246519;
    lapl[1]  = 0.12;

    //XC(mgga)(&mgga1, 1, rho, sigma, lapl, tau, 
    //	     &zk,  vrho, vsigma, vlapl, vtau, 
    //	     v2rho2, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau,
    //	     v2sigmalapl, v2sigmatau, v2lapltau);
    //XC(mgga)(&mgga2, 1, rho, sigma, lapl, tau, 
    //	     &zk2,  vrho2, vsigma2, vlapl2, vtau2, 
    // 	     NULL, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau,
    //	     v2sigmalapl, v2sigmatau, v2lapltau);
    //XC(mgga_exc)(&mgga2, 1, rho, sigma, lapl, tau, 
    //		 &zk2);
    XC(gga_exc_vxc)(&mgga1, 1, rho, sigma, 
		     &zk,  vrho, vsigma);
    XC(mgga_exc_vxc)(&mgga2, 1, rho, sigma, lapl, tau, 
		     &zk2,  vrho2, vsigma2, vlapl2, vtau2);

    fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", sigma[2], zk2*(rho[0]+rho[1]), vsigma2[2]);
  }

  XC(func_end)(&mgga1);
  XC(func_end)(&mgga2);
}

void test_neg_rho()
{
  xc_func_type func;
  double rho[5][2] = { {9.03897273e-06,-1.00463992e-06}, 
		       {8.48383564e-06,-3.51231267e-07}, 
		       {1.45740621e-08,-2.94546705e-09}, 
		       {2.62778445e-07, -1.00191745e-07}, 
		       {2.55745103e-06, -1.54789964e-06} };
  double sigma[5][3] = { {1.20122271e-08,4.83240746e-09,6.24774836e-09}, 
			 {1.54146602e-07,1.41584609e-07,1.36663204e-07}, 
			 {2.75312438e-08,2.75224049e-08,2.75135719e-08}, 
			 {1.90251649e-07,1.91241798e-07,1.92240989e-07}, 
			 {9.29562712e-09,7.83940082e-09, 8.05714636e-09} };
  double vsigma[5][3];
  double zk[5], vk[5][2];
  int i, func_id;

  for(func_id=1; func_id<1000; func_id++){
    if(xc_func_init(&func, func_id, XC_POLARIZED) != 0) continue;
    if(func_id == XC_LDA_C_2D_PRM || func_id == XC_GGA_X_LB) goto end;

    printf("\n%s:\n", func.info->name);

    switch(func.info->family){
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&func, 5, &rho[0][0], zk, &vk[0][0]);
      break;
    case XC_FAMILY_GGA:
      xc_gga_exc_vxc(&func, 5, &rho[0][0], &sigma[0][0], zk, &vk[0][0], &vsigma[0][0]);
      break;
    }

    switch(func.info->family){
    case XC_FAMILY_LDA:
      for(i=0; i<5; i+=1)
	printf("%.8e %.8e %.8e %.8e %.8e\n", 
	       rho[i][0], rho[i][1], zk[i], vk[i][0], vk[i][1]);
      break;
    case XC_FAMILY_GGA:
      for(i=0; i<5; i+=1)
	printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", 
	       rho[i][0], rho[i][1], sigma[i][0], sigma[i][1], sigma[i][2], 
	       zk[i], vk[i][0], vk[i][1], vsigma[i][0], vsigma[i][1], vsigma[i][2]);
      break;
    }

  end:
    XC(func_end)(&func);
  }
}


int main()
{
  //test_expi();
  //test_integration();
  //test_neg_rho();

  //test_lda();
  //test_enhance();
  test_gga(); 
  //test_ak13();
  //test_mgga();

  //printf("number = '%d'; key = '%s'", 25, XC(functional_get_name)(25));

  return 0;
}
