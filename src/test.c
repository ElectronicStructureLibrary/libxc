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
#include "util.h"

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
    result = integrate(func, NULL, a, b);
    printf("%lf %lf\n", b, result);
  }
}

void test_lda()
{
  XC(func_type) l1, l2, l3;
  int i;
  
  XC(func_init)(&l1, XC_LDA_C_GOMBAS, XC_UNPOLARIZED);
  XC(func_init)(&l2, XC_LDA_C_PW, XC_POLARIZED);
  XC(func_init)(&l3, XC_LDA_X, XC_UNPOLARIZED);

  //XC(lda_x_1d_set_params)(&l1, 0, 1.0);
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

    rho[0] = 0.01 + i/1000.0;
    rho[1] = 0.0;

    //rho[0] = 1.0/(2.0*rs);
    //rho[1] = 0.0;

    //dens = rho[0] + rho[1];

    XC(lda)(&l1, 1, rho, &ec1, vc1, fxc1, kxc1);
    XC(lda)(&l2, 1, rho, &ec2, vc2, NULL, NULL);
    //XC(lda_fxc_fd)(&l2, rho, fxc2);
    //XC(lda_kxc_fd)(&l2, rho, kxc2);

    //rho[0] = dens; rho[1] = 0.0;
    //XC(lda)(&l3, rho, &ec3, vc3, fxc3, kxc3);

    // printf("%e\t%e\t%e\n", dens, (fxc1[0]+2.0*fxc1[1]+fxc1[2])/4.0, fxc3[0]);
    // printf("%e\t%e\t%e\n", dens, (kxc1[0]+3.0*kxc1[1]+3.0*kxc1[2]+kxc1[3])/8.0, kxc3[0]);

    printf("%e\t%e\t%e\n", rho[0], fxc1[0], kxc1[0]);
  }
}

void test_gga()
{
  XC(func_type) gga, gga2;
  int i;

  //XC(func_init)(&gga, XC_GGA_C_AM05, XC_POLARIZED);
  XC(func_init)(&gga, XC_GGA_C_WL, XC_POLARIZED);
  //XC(func_init)(&gga2, XC_GGA_C_AM052, XC_UNPOLARIZED);
  
  for(i=0; i<=10000; i++){
    double rho[2], sigma[3], tau[2], lapl[2];
    double zk,   vrho[2],  vsigma[3];
    double zk2, vrho2[2], vsigma2[3];
    double v2rho2[3], v2sigma2[6], v2rhosigma[6];

    rho[0]   = 0.1 + i/10000.0;
    rho[1]   = 0.5;
    sigma[0] = 0.1;
    sigma[1] = 0.11;
    sigma[2] = 0.7;

    XC(gga)(&gga,  1, rho, sigma, &zk, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
    // XC(gga)(&gga2, 1, rho, sigma, &zk2, vrho2, vsigma2, v2rho2, v2rhosigma, v2sigma2);
    
    //fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", sigma[0], (rho[0])*zk, vsigma[0]);
    fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", rho[0], vrho[0], v2rho2[0]);
    //fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", rho[1], vrho[1], v2rho2[2]);
  }
}


//void 
//XC(mgga_c		     FLOAT *e, FLOAT *dedd, FLOAT *vsigma, FLOAT *dedtau);

void test_tpss()
{
  XC(func_type) tpss, tpss2;
  XC(func_type) agga;
  int i;

  XC(func_init)(&tpss, XC_MGGA_C_TPSS, XC_UNPOLARIZED);
  //XC(func_init)(&tpss2, XC_MGGA_C_TPSS_LARA, XC_UNPOLARIZED);
  //XC(mgga_c_tpss_init)(tpss2.mgga);
  
  for(i=0; i<=1000; i++){
    double rho[2], sigma[3], tau[2], lapl[2];
    double zk,   vrho[2],  vsigma[3],  vtau[2],  vlapl[2];
    double zk2, vrho2[2], vsigma2[3], vtau2[2], vlapl2[2];
    double v2rho2[3], v2sigma2[6], v2lapl2[3], v2tau2[3];
    double v2rhosigma[6], v2rholapl[3], v2rhotau[3];
    double v2sigmalapl[6], v2sigmatau[6], v2lapltau[3];

    rho[0]   = 0.1 + i/1000.0;
    rho[1]   = 0.15;
    sigma[0] = 0.1;
    sigma[1] = 0.11;
    sigma[2] = 0.7;
    tau[0]   = 0.03;
    tau[1]   = 0.15;
    lapl[0]  = 0.2;
    lapl[1]  = 0.12;

    XC(mgga)(&tpss, 1, rho,  sigma, lapl, tau, 
    	     &zk,  vrho, vsigma, vlapl, vtau, 
    	     NULL, v2sigma2, v2lapl2, v2tau2, v2rhosigma, v2rholapl, v2rhotau,
	     v2sigmalapl, v2sigmatau, v2lapltau);
    //XC(mgga_c_tpss_lara)(tpss2.mgga, rho, sigma, tau,
    //	    &zk2,  vrho2, vsigma2, vtau2);

    //XC(gga)(&agga, rho,  sigma,
    //&zk,  vrho, vsigma,
    //	    v2rho2, v2rhosigma, v2sigma2);
    fprintf(stderr, "%16.10lf\t%16.10lf\t%16.10lf\n", rho[0], zk, zk2);
  }
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
    xc_func_end(&func);
  }
}


int main()
{
  //test_integration();
  //test_neg_rho();

  //test_lda();
  test_gga();
  //test_tpss();

  return 0;
}
