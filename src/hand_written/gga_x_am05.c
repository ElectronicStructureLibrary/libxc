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

void XC(gga_x_am05_enhance)
  (const XC(func_type) *p, int order, FLOAT x, 
   FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  FLOAT ss, ss2, lam_x, dlam_x, d2lam_x, d3lam_x;
  FLOAT aux1, aux2, aux12, aux22;
  FLOAT ww, ww13, z_t, z_t2, z_tt, z_tt_aux, fx_b, xx, flaa_1, flaa_2, flaa;
  FLOAT dww, dz_t, dz_tt, dfx_b, dxx, dflaa_1, dflaa_2, dflaa;
  FLOAT d2ww, d2z_t, d2z_tt, d2fx_b, d2xx, d2flaa_1, d2flaa_2, d2flaa;
  FLOAT d3ww, d3z_t, d3z_tt, d3fx_b, d3xx, d3flaa_2, d3flaa;

  if(x < p->info->min_grad){
    *f    = 1.0;
    return;
  }

  ss  = X2S*x;
  ss2 = ss*ss;

  lam_x  = ss*SQRT(ss)/(2.0*SQRT(6.0));
  ww     = XC(lambert_w)(lam_x);
  ww13   = CBRT(ww);

  z_t    = (M_CBRT9/M_CBRT4)*ww13*ww13;
  z_t2   = z_t*z_t;

  /* This is equal to sqrt(t_zeta) * tt_zeta of the JCP*/
  z_tt_aux = d + z_t2;
  z_tt     = z_t * SQRT(SQRT(z_tt_aux));

  /* note that there is a factor of 2 missing in the JCP */
  fx_b   = M_PI/3.0*ss/z_tt;

  xx     = 1.0/(1.0 + alpha*ss2);

  flaa_1 = c*ss2 + 1.0;
  flaa_2 = c*ss2/fx_b + 1.0;
  flaa   = flaa_1/flaa_2;

  *f     = xx + (1.0 - xx)*flaa;

  if(order < 1) return;

  dlam_x  = 1.5*lam_x/ss;
  aux1    = 1.0 + ww;
  aux2    = lam_x*aux1;

  dww     = ww*dlam_x/aux2;
  dz_t    = M_CBRT2*dww/(M_CBRT3*ww13);
  dz_tt   = (2.0*d + 3.0*z_t2)*z_tt/(2.0*z_t*z_tt_aux);
  dfx_b   = M_PI/3.0*(z_tt - ss*dz_tt*dz_t)/(z_tt*z_tt);

  dxx     = -2.0*alpha*ss * xx*xx;
  dflaa_1 = 2.0*c*ss;
  dflaa_2 = DFRACTION(c*ss2, dflaa_1, fx_b, dfx_b);
  dflaa   = DFRACTION(flaa_1, dflaa_1, flaa_2, dflaa_2);

  *dfdx  = dxx*(1.0 - flaa) + dflaa*(1.0 - xx);
  *dfdx *= X2S;

  if(order < 2) return;

  aux12   = aux1*aux1;
  aux22   = aux2*aux2;

  d2lam_x  = 0.5*dlam_x/ss;
  d2ww     = ww*(-ww*(2.0 + ww)*dlam_x*dlam_x + aux12*lam_x*d2lam_x)/(aux22*aux1);
  d2z_t    = -M_CBRT2*(dww*dww - 3.0*ww*d2ww)/(3.0*M_CBRT3*ww*ww13);

  d2z_tt   = 3.0*z_t*(2.0*d + z_t2)*z_tt/(4.0*z_t*z_tt_aux*z_tt_aux);
  d2fx_b   = M_PI/3.0*(2.0*ss*dz_tt*dz_tt*dz_t*dz_t - z_tt*(dz_tt*(2.0*dz_t + ss*d2z_t) + ss*dz_t*dz_t*d2z_tt))/(z_tt*z_tt*z_tt);

  d2xx     = 2.0*alpha*(3.0*alpha*ss2 - 1.0) * xx*xx*xx;
  d2flaa_1 = 2.0*c;
  d2flaa_2 = D2FRACTION(c*ss2, dflaa_1, d2flaa_1, fx_b, dfx_b, d2fx_b);
  d2flaa   = D2FRACTION(flaa_1, dflaa_1, d2flaa_1, flaa_2, dflaa_2, d2flaa_2);

  *d2fdx2  = d2xx*(1.0 - flaa) - 2.0*dxx*dflaa + (1.0 - xx)*d2flaa;
  *d2fdx2 *= X2S*X2S;

  if(order < 3) return;

  d3lam_x  = -0.5*d2lam_x/ss;
  d3ww     = ww*(ww*dlam_x*(ww*(9.0 + 2.0*ww*(4.0 + ww))*dlam_x*dlam_x - 3.0*lam_x*aux12*(2.0 + ww)*d2lam_x) + 
		 lam_x*lam_x*aux12*aux12*d3lam_x)/(aux22*aux2*aux12);
  d3z_t    = M_CBRT2*(4.0*dww*dww*dww - 9.0*ww*dww*d2ww + 9.0*ww*ww*d3ww)/(9.0*M_CBRT3*ww*ww*ww13);

  d3z_tt   = -3.0*(-4.0*d*d + 4.0*d*z_t2 + z_t2*z_t2)*z_tt/(8.0*z_t*z_tt_aux*z_tt_aux*z_tt_aux);
  d3fx_b   = M_PI/3.0*(-6.0*ss*dz_t*dz_t*dz_t*dz_tt*dz_tt*dz_tt
		       +6.0*z_tt*dz_t*dz_tt*(dz_tt*(dz_t + ss*d2z_t) + ss*dz_t*dz_t*d2z_tt)
		       -z_tt*z_tt*(3.0*dz_t*(dz_t + ss*d2z_t)*d2z_tt + dz_tt*(3.0*d2z_t + ss*d3z_t) + ss*dz_t*dz_t*dz_t*d3z_tt))
    /(z_tt*z_tt*z_tt*z_tt);

  d3xx     = -24.0*alpha*alpha*ss*(alpha*ss2 - 1.0) * xx*xx*xx*xx;
  
  d3flaa_2 = D3FRACTION(c*ss2, dflaa_1, d2flaa_1, 0.0, fx_b, dfx_b, d2fx_b, d3fx_b);
  d3flaa   = D3FRACTION(flaa_1, dflaa_1, d2flaa_1, 0.0, flaa_2, dflaa_2, d2flaa_2, d3flaa_2);

  *d3fdx3  = d3xx*(1.0 - flaa) - 3.0*d2xx*dflaa + -3.0*dxx*d2flaa + (1.0 - xx)*d3flaa;
  *d3fdx3 *= X2S*X2S*X2S;

}

#define func XC(gga_x_am05_enhance)
