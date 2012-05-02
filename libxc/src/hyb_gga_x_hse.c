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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

#define XC_HYB_GGA_X_HSE03 427 /* the 2003 version of the screened hybrid HSE */
#define XC_HYB_GGA_X_HSE06 428 /* the 2006 version of the screened hybrid HSE */

typedef struct{
  FLOAT omega;
} hyb_gga_x_hse_params;


static void
hyb_gga_x_hse_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  assert(p->params == NULL);
  p->params = malloc(sizeof(hyb_gga_x_hse_params));

  /* value of omega in HSE */
  switch(p->info->number){
  case XC_HYB_GGA_X_HSE03:
    XC(hyb_gga_x_hse_set_params_)(p, 0.3);
    break;
  case XC_HYB_GGA_X_HSE06:
    XC(hyb_gga_x_hse_set_params_)(p, 0.2);
    break;
  default:
    fprintf(stderr, "Internal error in hyb_gga_x_hse\n");
    exit(1);
  }

  p->exx_coef = 0.25;
}


void 
XC(hyb_gga_x_hse_set_params)(XC(func_type) *p, FLOAT omega)
{
  assert(p != NULL && p->gga != NULL);
  XC(hyb_gga_x_hse_set_params_)(p->gga, omega);
}


void 
XC(hyb_gga_x_hse_set_params_)(XC(gga_type) *p, FLOAT omega)
{
  hyb_gga_x_hse_params *params;

  assert(p->params != NULL);
  params = (hyb_gga_x_hse_params *) (p->params);

  params->omega = omega;
}


#define HEADER 3

/* This implementation follows the one from espresso, that, in turn,
   follows the one of vasp. Analytic derivatives are only implemented
   in espresso though. These implementations can be found in:

   vasp: xclib_grad.F, MODULE wpbe, and in particular SUBROUTINE EXCHWPBE_R
   espresso: flib/functionals.f90, SUBROUTINE wpbe_analy_erfc_approx_grad

   very important details can be found in references:

   *) J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
   *) M Ernzerhof and JP Perdew, J. Chem. Phys. 109, 3313 (1998)
   *) J Heyd and GE Scuseria, J. Chem. Phys. 120, 7274 (2004)
*/

static inline void 
func(const XC(gga_type) *p, int order, FLOAT x, FLOAT ds,
     FLOAT *f, FLOAT *dfdx, FLOAT *lvrho)
{
  static const FLOAT AA=1.0161144, BB=-0.37170836, CC=-0.077215461, DD=0.57786348, EE=-0.051955731;
  static const FLOAT m89=-8.0/9.0;

  /* Cutoff criterion below which to use polynomial expansion */
  static const FLOAT EGscut=0.08, wcutoff=14, expfcutoff=700.0;

  /* parameters for the re-scaling of s */
  static const FLOAT strans=8.3, smax=8.5728844, sconst=18.79622316;

  FLOAT omega, kF, ww, ww2, ww3, ww4, ww5, ww6, ww7, ww8;
  FLOAT ss, ss2, ss3, ss4, ss5, ss6;
  FLOAT DHs, DHs2, DHs3, DHs4, DHs72, DHs92;
  FLOAT f94Hs2_A, DHsw, DHsw2, DHsw52, DHsw72;
  FLOAT Hsbw, Hsbw2, Hsbw3, Hsbw4, Hsbw12, Hsbw32, Hsbw52, Hsbw72;
  FLOAT DHsbw, DHsbw2, DHsbw3, DHsbw4, DHsbw5, DHsbw12, DHsbw32, DHsbw52, DHsbw72, DHsbw92;
  FLOAT HsbwA94, HsbwA942, HsbwA943, HsbwA945, HsbwA9412;
  FLOAT eb1, H, F, EG;
  FLOAT term1, term2, term3, term4, term5, t10, piexperf, expei;

  /* x   = 10.3728/X2S; */
  /* ds  = 0.1234; */

  assert(p->params != NULL);
  omega = ((hyb_gga_x_hse_params *)(p->params))->omega;


  kF  = POW(3.0*M_PI*M_PI*ds, 1.0/3.0);
  ww  = omega/kF;
  ww2 = ww*ww; ww3 = ww*ww2; ww4 = ww*ww3; ww5 = ww*ww4; ww6 = ww*ww5; ww7 = ww*ww6; ww8 = ww*ww7;
  /* printf("kF = %le; w = %le\n", kF, ww); */

  /*  Rescaling the s values to ensure the Lieb-Oxford bound for s>8.3 */
  ss  = (X2S*x < strans) ? X2S*x : smax - sconst/(X2S*X2S*x*x);
  ss2 = ss*ss; 
  ss3 = ss*ss2; 
  ss4 = ss*ss3; 
  ss5 = ss*ss4; 
  ss6 = ss*ss5;

  /* first let us calculate H(s) */
  {
    static const FLOAT Ha1=0.00979681, Ha2=0.0410834, Ha3=0.187440, Ha4=0.00120824, Ha5=0.0347188;
    FLOAT Hnum, Hden;

    Hnum = Ha1*ss2 + Ha2*ss4;
    Hden = 1.0 + Ha3*ss4 + Ha4*ss5 + Ha5*ss6;

    H = Hnum/Hden;
  }
  /* printf("H = %le\n", H); */

  /* useful variables for what comes next */
  DHs   = DD + H*ss2; 
  DHs2  = DHs*DHs; 
  DHs3  = DHs2*DHs; 
  DHs4  = DHs3*DHs;
  DHs72 = DHs3*SQRT(DHs); 
  DHs92 = DHs72*DHs;

  f94Hs2_A = 9.0*H*ss2/(4.0*AA);

  DHsw   = DHs + ww2;
  DHsw2  = DHsw*DHsw; 
  DHsw52 = SQRT(DHsw)*DHsw2; 
  DHsw72 = DHsw52*DHsw;

  eb1 = (ww < wcutoff) ? 1.455915450052607 : 2.0;

  Hsbw   = ss2*H + eb1*ww2; 
  Hsbw2  = Hsbw*Hsbw; 
  Hsbw3  = Hsbw2*Hsbw; 
  Hsbw4  = Hsbw3*Hsbw;
  Hsbw12 = SQRT(Hsbw); 
  Hsbw32 = Hsbw12*Hsbw; 
  Hsbw52 = Hsbw32*Hsbw; 
  Hsbw72 = Hsbw52*Hsbw;

  DHsbw   = DD + ss2*H + eb1*ww2;
  DHsbw2  = DHsbw*DHsbw; 
  DHsbw3  = DHsbw2*DHsbw; 
  DHsbw4  = DHsbw3*DHsbw; 
  DHsbw5  = DHsbw4*DHsbw;
  DHsbw12 = SQRT(DHsbw); 
  DHsbw32 = DHsbw12*DHsbw; 
  DHsbw52 = DHsbw32*DHsbw; 
  DHsbw72 = DHsbw52*DHsbw;
  DHsbw92 = DHsbw72*DHsbw;

  HsbwA94   = 9.0*Hsbw/(4.0*AA);
  HsbwA942  = HsbwA94*HsbwA94;
  HsbwA943  = HsbwA942*HsbwA94;
  HsbwA945  = HsbwA943*HsbwA942;
  HsbwA9412 = SQRT(HsbwA94);

  /* now we calculate F(s) */
  {
    static const FLOAT Fc1=4.0*AA*AA/(9.0*CC) + (BB - AA*DD)/CC, Fc2=-4.0/(3.0*36.0*CC);

    F = Fc1*H + Fc2;
  }
  /* printf("F = %le\n", F); */

  /* and now G(s) */
  if(ss > EGscut){
    FLOAT Ga, Gb;

    Ga = M_SQRTPI*(15.0*EE + 6.0*CC*(1.0 + F*ss2)*DHs + 4.0*BB*DHs2 + 8.0*AA*DHs3)/(16.0*DHs72)
      - (3.0*M_PI/4.0)*sqrt(AA)*exp(f94Hs2_A)*(1.0 - erf(sqrt(f94Hs2_A)));
    Gb = 15.0*M_SQRTPI*ss2/(16.0*DHs72);

    EG = -(3.0*M_PI/4.0 + Ga)/Gb;
  }else{
    static const FLOAT EGa1=-0.02628417880, EGa2=-0.07117647788, EGa3=0.08534541323;

    EG = EGa1 + EGa2*ss2 + EGa3*ss4;
  }
  /* printf("EG = %le\n", EG); */

  /* Calculate the terms needed in any case */
  term2 = (DHs2*BB + DHs*CC + 2.0*EE + DHs*ss2*CC*F + 2.0*ss2*EG)/(2.0*DHs3);
  term3 = -ww*(4.0*DHsw2*BB + 6.0*DHsw*CC + 15.0*EE + 6.0*DHsw*ss2*CC*F + 15.0*ss2*EG)/(8.0*DHs*DHsw52);
  term4 = -ww3*(DHsw*CC + 5.0*EE + DHsw*ss2*CC*F + 5.0*ss2*EG)/(2.0*DHs2*DHsw52);
  term5 = -ww5*(EE + ss2*EG)/(DHs3*DHsw52);
  
  /* printf("%le %le %le %le\n", term2, term3, term4, term5); */

  if((ss > 0.0) || (ww > 0.0)){
    t10 = AA*LOG(Hsbw / DHsbw)/2.0;
  }

  /* Calculate exp(x)*f(x) depending on size of x */
  if(HsbwA94 < expfcutoff){
    piexperf = M_PI*exp(HsbwA94)*erfc(HsbwA9412);
    expei    = exp(HsbwA94)*(-expint(HsbwA94));
  }else{
    static const FLOAT expei1=4.03640, expei2=1.15198, expei3=5.03627, expei4=4.19160;

    piexperf = M_PI*(1.0/(M_SQRTPI*HsbwA9412) - 1.0/(2.0*SQRT(M_PI*HsbwA943))+ 3.0/(4.0*SQRT(M_PI*HsbwA945)));
    expei  = - (1.0/HsbwA94)*(HsbwA942 + expei1*HsbwA94 + expei2)/(HsbwA942 + expei3*HsbwA94 + expei4);
  }
  /* printf("exp = %le, %le\n", piexperf, expei); */

  if (ww == 0.0){ /* Fall back to original expression for the PBE hole */
    FLOAT t1 = -AA*expei/2.0;
    if(ss > MIN_GRAD){
      term1 = t1 + t10;
      *f = m89*(term1 + term2);
    }else{
      *f = 1.0;
    }

  }else if(ww > wcutoff){ /* Use simple gaussian approximation for large w */
    term1   = -AA*(expei + LOG(DHsbw) - LOG(Hsbw))/2.0;
    *f = m89*(term1 + term2 + term3 + term4 + term5);

  }else{ /*  For everything else use the full blown expression */

    static const FLOAT ea1=-1.128223946706117, ea2=1.452736265762971, ea3=-1.243162299390327,
      ea4=0.971824836115601, ea5=-0.568861079687373, ea6=0.246880514820192, ea7=-0.065032363850763,
      ea8=0.008401793031216;

    FLOAT AA2, AA3, AA4, AA12, AA32, AA52, AA72;
    FLOAT np1, np2, t1, f1, f2, f3, f4, f5, f6, f7, f8, f9, t2t9;

    AA2  = AA*AA;
    AA3  = AA2*AA;
    AA4  = AA3*AA;
    AA12 = SQRT(AA);
    AA32 = AA12*AA;
    AA52 = AA32*AA;
    AA72 = AA52*AA;

    np1 = -1.5*ea1*AA12*ww + 27.0*ea3*ww3/(8.0*AA12) - 243.0*ea5*ww5/(32.0*AA32) + 2187.0*ea7*ww7/(128.0*AA52);
    np2 = -AA + 9.0*ea2*ww2/4.0 - 81.0*ea4*ww4/(16.0*AA);

    t1 = 0.5*(np1*piexperf + np2*expei);

    f2 = 0.5*ea1*M_SQRTPI*AA/DHsbw12;
    f3 = 0.5*ea2*AA/DHsbw;
    f4 = ea3*M_SQRTPI*(-9.0/(8.0*Hsbw12) + 0.25*AA/DHsbw32);
    f5 = (ea4/128.0)*(-144.0/Hsbw + 64.0*AA/DHsbw2);
    f6 = ea5*(3.0*M_SQRTPI*(3.0*DHsbw52*(9.0*Hsbw - 2.0*AA)
			    + 4.0*Hsbw32*AA2))/(32.0*DHsbw52*Hsbw32*AA);
    f7 = ea6*((32.0*AA/DHsbw3 + (-36.0 + 81.0*ss2*H/AA)/Hsbw2))/32.0;
    f8 = ea7*(-3.0*M_SQRTPI*(-40.0*Hsbw52*AA3 + 9.0*DHsbw72*(27.0*Hsbw2 - 6.0*Hsbw*AA + 4.0*AA2)))/(128.0*DHsbw72*Hsbw52*AA2);
    f9 = (324.0*ea6*eb1*DHsbw4*Hsbw*AA + ea8*(384.0*Hsbw3*AA3 + DHsbw4*(-729.0*Hsbw2 + 324.0*Hsbw*AA - 288.0*AA2)))/(128.0*DHsbw4*Hsbw3*AA2);

    /*
    printf("t1 = %le\n", t1);
    printf("f2 = %le\n", f2);
    printf("f3 = %le\n", f3);
    printf("f4 = %le\n", f4);
    printf("f5 = %le\n", f5);
    printf("f6 = %le\n", f6);
    printf("f7 = %le\n", f7);
    printf("f8 = %le\n", f8);
    printf("f9 = %le\n", f9);
    */

    t2t9  = f2*ww + f3*ww2 + f4*ww3 + f5*ww4 + f6*ww5 + f7*ww6 + f8*ww7 + f9*ww8;

    term1 = t1 + t2t9 + t10;

    *f = m89*(term1 + term2 + term3 + term4 + term5);

    /* printf("F1 = %le %le %le\n", term1, term5, m89); */
  }
  /* printf("F = %le\n", *f); */
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_hyb_gga_x_hse03) = {
  XC_HYB_GGA_X_HSE03,
  XC_EXCHANGE,
  "HSE03",
  XC_FAMILY_HYB_GGA,
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_x_hse_init,
  NULL, NULL, 
  work_gga_x
};

const XC(func_info_type) XC(func_info_hyb_gga_x_hse06) = {
  XC_HYB_GGA_X_HSE06,
  XC_EXCHANGE,
  "HSE06",
  XC_FAMILY_HYB_GGA,
  "J Heyd, GE Scuseria, and M Ernzerhof, J. Chem. Phys. 118, 8207 (2003)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC,
  1e-32, 1e-32, 0.0, 1e-32,
  hyb_gga_x_hse_init,
  NULL, NULL, 
  work_gga_x
};
