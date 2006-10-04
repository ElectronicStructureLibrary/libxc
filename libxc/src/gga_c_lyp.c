#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

void gga_c_lyp(void *p_, double *rho_, double *sigma_,
	       double *e, double *vrho, double *vsigma)
{
  gga_type *p = p_;
  
  static double AA = 0.04918, BB = 0.132, cc = 0.2533, 
    dd = 0.349, ee = 36.462398978764767321; /* ee = 8*2^(2/3)*e */

  double rho[2], sigma[2], rhot, sigmat;
  double sfact, rhot13, rhot43, rho83[2], ZZ, delta, omega;
  double dZZdr, ddeltadr, domegadr;
  double t1, t2, t3, t4, t5, t6;
  int is;

  /* convert input to spin-polarized case */
  if(p->nspin == XC_POLARIZED){
    sfact    = 1.0;
    rho[0]   = rho_[0];
    rho[1]   = rho_[1];
    sigma[0] = sigma_[0];
    sigma[1] = sigma_[2];
    rhot     = rho_[0] + rho_[1];
    sigmat   = sigma_[0] + 2.0*sigma_[1] + sigma_[2];
  }else{
    sfact    = 2.0;
    rho[0]   = rho_[0]/2.0;
    rho[1]   = rho_[0]/2.0;
    sigma[0] = sigma_[0]/4.0;
    sigma[1] = sigma_[0]/4.0;
    rhot     = rho_[0];
    sigmat   = sigma_[0];    
  }

  /* some handy functions of the total density */
  rhot13   = pow(rhot, 1.0/3.0);
  rhot43   = rhot*rhot13;
  rho83[0] = pow(rho[0], 8.0/3.0);
  rho83[1] = pow(rho[1], 8.0/3.0);

  ZZ     = rhot13/(rhot13 + dd);
  delta  = (cc + dd*ZZ)/rhot13;
  omega  = exp(-cc/rhot13) * ZZ * pow(rhot, -11.0/3.0);

  /* and their derivatives */
  dZZdr    = dd*ZZ*ZZ/(3.0*rhot43);
  ddeltadr = -(cc + dd*ZZ - dd*dd*ZZ*ZZ/rhot13)/(3.0*rhot43);
  domegadr = omega*ZZ*(cc*dd + (cc-10.0*dd)*rhot13 - 11.0*rhot13*rhot13)/(3.0*rhot43*rhot13);

  /* we now add the terms of the xc energy one by one */

  /* t1 */
  {
    double aux1;

    aux1 = -4.0*AA/rhot;

    t1 = aux1*rho[0]*rho[1]*ZZ;

    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;

      vrho[is] = aux1*(rho[js]*ZZ - rho[0]*rho[1]*(ZZ/rhot - dZZdr));
    }
  }

  
  /* t2 */
  {
    double aux1, aux2, aux3, aux4;

    aux1 = -AA*BB;
    aux2 = rho[0]*rho[1];
    aux3 = (47.0 - 7.0*delta)/18.0;
    aux4 = 2.0/3.0*rhot*rhot;

    t2 = aux1*omega*sigmat*(aux2*aux3 - aux4);

    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;

      vrho[is] += aux1*
	(domegadr*sigmat*(aux2*aux3 - aux4) +
	 omega*sigmat*(rho[js]*aux3 - aux2*7.0/18.0*ddeltadr - 4.0/3.0*rhot));
    }

    {
      double aux5 = aux1*omega*(aux2*aux3 - aux4);

      vsigma[0] += aux5/sfact;
      if(p->nspin == XC_POLARIZED){
	vsigma[1] += 2.0*aux5;
	vsigma[2] +=     aux5;
      }else{
	vsigma[0] += aux5/2.0;
      }
    }
  }


  /* t3 */
  {
    double aux1, aux2;
    aux1 = -AA*BB;
    aux2 = ee*(rho83[0] + rho83[1]);

    t3 = aux1*omega*rho[0]*rho[1]*aux2;
    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;

      vrho[is] += aux1*domegadr*rho[0]*rho[1]*aux2;
      vrho[is] += aux1* omega * ee
	* rho[js]*(11.0/3.0*rho83[is] + rho83[js]);
    }
  }


  /* t4 */
  {
    double aux1, aux2, aux3, aux4;

    aux1 = AA*BB;
    aux2 = rho[0]*rho[1];
    aux3 = 5.0/2.0 - delta/18.0;
    aux4 = sigma[0] + sigma[1];

    t4 = aux1*omega*aux2*aux3*aux4;

    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;
      
      vrho[is] += aux1*aux4*
	(domegadr*aux2*aux3 + omega*rho[js]*aux3 -
	 omega*aux2*ddeltadr/18.0);
    }

    {
      double aux5 = aux1*omega*aux2*aux3;

      vsigma[0] += aux5/sfact;
      if(p->nspin == XC_POLARIZED)
	vsigma[2] += aux5;
    }
  }


  /* t5 */
  {
    double aux1, aux2, aux3, aux4;

    aux1 = AA*BB;
    aux2 = rho[0]*rho[1]/(9.0*rhot);
    aux3 = delta - 11.0;
    aux4 = rho[0]*sigma[0] + rho[1]*sigma[1];

    t5 = aux1*omega*aux2*aux3*aux4;

    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;

      vrho[is] += aux1*
	((domegadr*aux2 + omega*(rho[js]*rhot - rho[0]*rho[1])/(9.0*rhot*rhot))*aux3*aux4 +
	 omega*aux2*(ddeltadr*aux4 + aux3*sigma[is]));

      vsigma[is==0 ? 0 : 2] += aux1*omega*aux2*aux3*rho[is]/sfact;
    }
  }

  
  /* t6 */
  {
    double aux1, aux2, aux3;

    aux1 = -AA*BB;
    aux2 = 2.0/3.0*rhot*rhot*(sigma[0] + sigma[1]);
    aux3 = rho[0]*rho[0]*sigma[1] + rho[1]*rho[1]*sigma[0];

    t6 = aux1*omega*(aux2 - aux3);
    for(is=0; is<p->nspin; is++){
      int js = (is==0) ? 1 : 0;

      vrho[is] += aux1*
	(domegadr*(aux2 - aux3) + 
	 omega*(4.0/3.0*rhot*(sigma[0] + sigma[1]) - 2*rho[is]*sigma[js]));

      vsigma[is==0 ? 0 : 2] += aux1*omega*(2.0/3.0*rhot*rhot - rho[js]*rho[js])/sfact;
    }
  }

  /* we add all contributions to the total energy */
  *e = (t1 + t2 + t3 + t4 + t5 + t6)/rhot;
}


const func_type func_gga_c_lyp = {
  XC_GGA_C_LYP,
  XC_CORRELATION,
  "Lee, Yang & Parr",
  XC_FAMILY_GGA,
  "C. Lee, W. Yang and R.G. Parr, Phys. Rev. B 37, 785 (1988)\n"
  "B. Miehlich, A. Savin, H. Stoll and H. Preuss, Chem. Phys. Letters 157, 200 (1989)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  gga_c_lyp
};
