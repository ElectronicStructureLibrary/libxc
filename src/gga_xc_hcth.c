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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "util.h"

#define XC_GGA_XC_HCTH_93  161 /* HCTH functional fitted to  93 molecules  */
#define XC_GGA_XC_HCTH_120 162 /* HCTH functional fitted to 120 molecules  */
#define XC_GGA_XC_HCTH_147 163 /* HCTH functional fitted to 147 molecules  */
#define XC_GGA_XC_HCTH_407 164 /* HCTH functional fitted to 147 molecules  */

static void gga_xc_hcth_init(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  p->lda_aux = (xc_lda_type *) malloc(sizeof(xc_lda_type));
  xc_lda_init(p->lda_aux, XC_LDA_C_PW, XC_POLARIZED);
}

static void gga_xc_hcth_end(void *p_)
{
  xc_gga_type *p = (xc_gga_type *)p_;

  free(p->lda_aux);
}

void func_g(int func, int type, double s, double *g, double *dg, double *ldg)
{
  const double c[4][3][5] = {
    {      /* HCTH/93 */
      {1.09320,  -0.744056,    5.59920,   -6.78549,   4.49357}, /* X   */
      {0.222601, -0.0338622,  -0.0125170, -0.802496,  1.55396}, /* Css */
      {0.729974,  3.35287,   -11.5430,     8.08564,  -4.47857}  /* Cab */
    }, {   /* HCTH/120 */
      {1.09163,  -0.747215,  5.07833,  -4.10746,   1.17173},    /* X   */
      {0.489508, -0.260699,  0.432917, -1.99247,   2.48531},    /* Css */
      {0.514730,  6.92982, -24.7073,   23.1098,  -11.3234 }     /* Cab */
    }, {   /* HCTH/147 */
      {1.09025, -0.799194,   5.57212, -5.86760,  3.04544 },     /* X   */
      {0.562576, 0.0171436, -1.30636,  1.05747,  0.885429},     /* Css */
      {0.542352, 7.01464,  -28.3822,  35.0329, -20.4284  },     /* Cab */
    }, {   /* HCTH/407 */
      {1.08184, -0.518339,  3.42562, -2.62901,  2.28855},       /* X   */
      {1.18777, -2.40292,   5.61741, -9.17923,  6.24798},       /* Css */
      {0.589076, 4.42374, -19.2218,  42.5721, -42.0052 }        /* Cab */
    }
  };
  const double gamma[3] = {
    0.004, 0.2, 0.006
  };

  double s2, dd, x, dx;
  const double *cc;

  s2 = s*s;
  dd = (1.0 + gamma[type]*s2);
  x  = gamma[type] * s2/dd;
  dx = gamma[type] * 2.0*s/(dd*dd);

  cc = c[func][type];

  *g   = cc[0] + x*(cc[1] + x*(cc[2] + x*(cc[3] + x*cc[4])));
  *dg  = cc[1] + x*(2.0*cc[2] + x*(3.0*cc[3] + x*4.0*cc[4]));
  *dg *= dx;
  *ldg = cc[1]*gamma[type]; /* times 2s */
}

static void 
gga_xc_hcth(void *p_, double *rho, double *sigma,
	    double *e, double *vrho, double *vsigma)
{
  xc_gga_type *p = p_;

  double dens, mrho[2], ecunif, vcunif[2], x_avg, x[2];
  double sfact;
  int func, is;

  switch(p->info->number){
  case XC_GGA_XC_HCTH_120: func = 1; break;
  case XC_GGA_XC_HCTH_147: func = 2; break;
  case XC_GGA_XC_HCTH_407: func = 3; break;
  default:                 func = 0; /* XC_GGA_XC_HCTH_93 */
  }

  *e = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact   = 1.0;
    mrho[0] = rho[0];
    mrho[1] = rho[1];
    dens    = rho[0] + rho[1];
    vsigma[1] = 0.0;
  }else{
    sfact   = 2.0;
    mrho[0] = rho[0]/2.0;
    mrho[1] = mrho[0];
    dens    = rho[0];
  }

  xc_lda_vxc(p->lda_aux, mrho, &ecunif, vcunif);
  ecunif *= dens;

  x_avg = 0.0;
  for(is=0; is<p->nspin; is++){
    double mrho2[2], gdm, ds, rho13;
    double g_x, dg_x, ldg_x, g_ss, dg_ss, ldg_ss, e_x, e_ss;
    double ecunif_s, vcunif_s[2];
    int js = is==0 ? 0 : 2;

    vrho[is]   = 0.0;
    vsigma[js] = 0.0;
    if(rho[is] < MIN_DENS) continue;

    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = pow(ds, 1.0/3.0);
    x[is] = gdm/(ds*rho13);
    x_avg+= sfact*0.5*x[is]*x[is];

    /* the exchange term */
    func_g(func, 0, x[is], &g_x, &dg_x, &ldg_x);

    e_x       = -sfact*X_FACTOR_C*(ds*rho13);
    vrho[is] += -4.0/3.0*X_FACTOR_C*rho13*(g_x - dg_x*x[is]);

    /* the ss term */
    mrho2[0] = mrho[is];
    mrho2[1] = 0.0;
    xc_lda_vxc(p->lda_aux, mrho2, &ecunif_s, vcunif_s);
    func_g(func, 1, x[is], &g_ss, &dg_ss, &ldg_ss);

    e_ss      = sfact*ds*ecunif_s;
    vrho[is] += vcunif_s[0]*g_ss - 4.0/3.0*ecunif_s*dg_ss*x[is];

    (*e) += e_x*g_x + e_ss*g_ss;

    if(gdm>MIN_GRAD)
      vsigma[js] = (e_x*dg_x + e_ss*dg_ss) * x[is]/(2.0*sigma[js]);
    else
      vsigma[js] = (e_x*ldg_x + e_ss*ldg_ss) / (sfact*sfact*ds*rho13);

    /* correct to get e_ab and v_ab */
    ecunif     -= sfact*ds*ecunif_s;
    vcunif[is] -= vcunif_s[0];    
  }

  { /* now the ab term */
    double g_ab, dg_ab, ldg_ab;

    x_avg = sqrt(x_avg);
    func_g(func, 2, x_avg, &g_ab, &dg_ab, &ldg_ab);
    (*e) += ecunif*g_ab;

    for(is=0; is<p->nspin; is++){
      double dd;
      int js = is==0 ? 0 : 2;

      vrho[is] += vcunif[is]*g_ab;

      dd = pow(dens, 4.0/3.0);
      if(x_avg*dd > MIN_GRAD*MIN_GRAD && mrho[is] > MIN_DENS){
	vrho[is]   -= 4.0/3.0*(ecunif/mrho[is])*dg_ab*x[is]*x[is]/(2.0*x_avg);
	vsigma[js] += ecunif*dg_ab*pow(mrho[is], -8.0/3.0)/(sfact*4.0*x_avg);
      }
    }
  }

  *e /= dens; /* we want energy per particle */
}


const xc_func_info_type func_info_gga_xc_hcth_93 = {
  XC_GGA_XC_HCTH_93,
  XC_EXCHANGE_CORRELATION,
  "HCTH/93",
  XC_FAMILY_GGA,
  "FA Hamprecht, AJ Cohen, DJ Tozer, and NC Handy, J. Chem. Phys. 109 6264 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_hcth_init, 
  gga_xc_hcth_end, 
  NULL,
  gga_xc_hcth
};

const xc_func_info_type func_info_gga_xc_hcth_120 = {
  XC_GGA_XC_HCTH_120,
  XC_EXCHANGE_CORRELATION,
  "HCTH/120",
  XC_FAMILY_GGA,
  "AD Boese, NL Doltsinis, NC Handy, and M Sprik, J. Chem. Phys. 112 1670 (2000)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_hcth_init, 
  gga_xc_hcth_end, 
  NULL,
  gga_xc_hcth
};

const xc_func_info_type func_info_gga_xc_hcth_147 = {
  XC_GGA_XC_HCTH_147,
  XC_EXCHANGE_CORRELATION,
  "HCTH/147",
  XC_FAMILY_GGA,
  "AD Boese, NL Doltsinis, NC Handy, and M Sprik, J. Chem. Phys. 112 1670 (2000)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_hcth_init, 
  gga_xc_hcth_end, 
  NULL,
  gga_xc_hcth
};

const xc_func_info_type func_info_gga_xc_hcth_407 = {
  XC_GGA_XC_HCTH_407,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407",
  XC_FAMILY_GGA,
  "AD Boese, and NC Handy, J. Chem. Phys. 114 5497 (2001)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  gga_xc_hcth_init, 
  gga_xc_hcth_end, 
  NULL,
  gga_xc_hcth
};
