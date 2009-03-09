/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <math.h>
#include <stdio.h>

#include "util.h"
#include "xc.h"
#include "util.h"

void test_lda()
{
  XC(lda_type) l1, l2, l3;
  int i;
  
  XC(lda_x_init)(&l1, XC_POLARIZED, 3, XC_NON_RELATIVISTIC);
  XC(lda_init)(&l2, XC_LDA_C_PZ, XC_POLARIZED);
  XC(lda_init)(&l3, XC_LDA_C_PZ, XC_UNPOLARIZED);

  for(i=0; i<1000; i++){
    double dens, rs, zeta, rho[2];
    double ec1, vc1[2], fxc1[3], kxc1[4];
    double ec2, vc2[2], fxc2[3], kxc2[4];
    double ec3, vc3[2], fxc3[3], kxc3[4];
    
    rs   = 1.0;
    zeta = -1.0 + 2.0*i/1000.0;

    dens = 1.0/(4.0/3.0*M_PI*POW(rs,3)); /* 3D */
    //dens = 1.0/(2.0*rs); /* 1D */

    //rho[0] = dens*(1.0 + zeta)/2.0;
    //rho[1] = dens*(1.0 - zeta)/2.0;

    rho[0] = 0.01 + i/1000.0;
    rho[1] = 0.21;

    XC(lda)(&l2, rho, &ec1, vc1, fxc1, kxc1);
    XC(lda_kxc_fd)(&l2, rho, kxc2);

    //rho[0] = dens; rho[1] = 0.0;
    //XC(lda)(&l3, rho, &ec3, vc3, fxc3, kxc3);

    // printf("%e\t%e\t%e\n", dens, (fxc1[0]+2.0*fxc1[1]+fxc1[2])/4.0, fxc3[0]);
    // printf("%e\t%e\t%e\n", dens, (kxc1[0]+3.0*kxc1[1]+3.0*kxc1[2]+kxc1[3])/8.0, kxc3[0]);

    printf("%e\t%e\t%e\n", rho[0], kxc1[3], kxc2[3]);
  }
}

void test_tpss()
{
  XC(mgga_type) tpss;
  XC(gga_type) agga;
  int i;

  //XC(mgga_init)(&tpss, XC_MGGA_X_LTA, XC_POLARIZED);
  XC(gga_init)(&agga, XC_GGA_XC_B97, XC_UNPOLARIZED);
  
  for(i=0; i<1000; i++){
    double rho[2], sigma[3], tau[2];
    double zk, vrho[2], vsigma[3], vtau[2];
    double v2rho2[3], v2rhosigma[6], v2sigma2[6], v2rhotau[4], v2tausigma[6], v2tau2[3];

    rho[0]   = 0.23;
    rho[1]   = 0.15;
    sigma[0] = 0.01 + i/1000.0;
    sigma[1] = 0.11;
    sigma[2] = 0.7;
    tau[0]   = 0.22;
    tau[1]   = 0.15;

    //XC(mgga)(&agga, rho,  sigma,  tau, 
    //	     &zk,  vrho, vsigma, vtau, 
    //	     v2rho2, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
    XC(gga)(&agga, rho,  sigma,
	    &zk,  vrho, vsigma,
	    v2rho2, v2rhosigma, v2sigma2);
    printf("%16.10lf\t%16.10lf\t%16.10lf\t%16.10lf\n", sigma[0], (rho[0]+rho[1])*zk, vrho[0], v2rhosigma[0]);
  }
}

int main()
{
  test_lda();

  return 0;
}
