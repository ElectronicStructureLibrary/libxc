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
  double rs;
  
  XC(lda_x_init)(&l1, XC_POLARIZED, 3, XC_NON_RELATIVISTIC);
  XC(lda_init)(&l2, XC_LDA_C_PZ, XC_POLARIZED);
  XC(lda_init)(&l3, XC_LDA_C_VWN, XC_POLARIZED);

  for(rs=10; rs>=0.1; rs-=0.01){
    double dens, zeta, rho[2];
    double ec1, vc1[2], fxc1[4];
    double ec2, vc2[2], fxc2[4];
    
    dens = 1.0/(4.0/3.0*M_PI*POW(rs,3));
    /* zeta = 1;
       rho[0] = dens*(1.0 + zeta)/2.0;
       rho[1] = dens*(1.0 - zeta)/2.0; */

    rho[0] = dens;
    rho[1] = 0;
    dens   = (rho[0] + rho[1]);
    zeta   = (rho[0] - rho[1])/dens;

    XC(lda_vxc)(&l2, rho, &ec1, vc1);
    XC(lda_vxc)(&l3, rho, &ec2, vc2);
    XC(lda_fxc)(&l1, rho, fxc1);
    XC(lda_fxc)(&l3, rho, fxc2);

    printf("%e\t%e\t%e\t%e\t%e\n", dens, fxc1[0], fxc1[1], fxc1[2], fxc1[3]);
    
  }
}

void test_tpss()
{
  XC(mgga_type) tpss;
  int i;

  XC(mgga_init)(&tpss, XC_MGGA_X_M06L, XC_POLARIZED);
  
  for(i=0; i<1000; i++){
    double rho[2], sigma[3], tau[2];
    double zk, vrho[2], vsigma[3], vtau[2];
    double v2rho2[3], v2rhosigma[6], v2sigma2[6], v2rhotau[4], v2tausigma[6], v2tau2[3];

    rho[0]   = 0.01 + i/1000.0;
    rho[1]   = 0.23;
    sigma[0] = 0.23;
    sigma[1] = 0.11;
    sigma[2] = 0.7;
    tau[0]   = 0.44;
    tau[1]   = 0.15;

    XC(mgga)(&tpss, rho,  sigma,  tau, 
	     &zk,  vrho, vsigma, vtau, 
	     NULL, v2rhosigma, v2sigma2, v2rhotau, v2tausigma, v2tau2);
    printf("%16.10lf\t%16.10lf\t%16.10lf\n", rho[0], (rho[0]+rho[1])*zk, vrho[0]);
  }
}

int main()
{
  test_tpss();

  return 0;
}
