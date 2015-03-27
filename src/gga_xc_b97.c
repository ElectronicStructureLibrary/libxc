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
#include <assert.h>
#include <stdlib.h>
#include "util.h"

#define XC_GGA_XC_HCTH_93  161 /* HCTH functional fitted to  93 molecules  */
#define XC_GGA_XC_HCTH_120 162 /* HCTH functional fitted to 120 molecules  */
#define XC_GGA_XC_HCTH_147 163 /* HCTH functional fitted to 147 molecules  */
#define XC_GGA_XC_HCTH_407 164 /* HCTH functional fitted to 407 molecules  */
#define XC_GGA_XC_B97      167 /* Becke 97                                 */
#define XC_GGA_XC_B97_1    168 /* Becke 97-1                               */
#define XC_GGA_XC_B97_2    169 /* Becke 97-2                               */
#define XC_GGA_XC_B97_D    170 /* Grimme functional to be used with C6 vdW term */
#define XC_GGA_XC_B97_K    171 /* Boese-Martin for Kinetics                */
#define XC_GGA_XC_B97_3    172 /* Becke 97-3                               */
#define XC_GGA_XC_SB98_1a  176 /* Schmider-Becke 98 parameterization 1a    */
#define XC_GGA_XC_SB98_1b  177 /* Schmider-Becke 98 parameterization 1b    */
#define XC_GGA_XC_SB98_1c  178 /* Schmider-Becke 98 parameterization 1c    */
#define XC_GGA_XC_SB98_2a  179 /* Schmider-Becke 98 parameterization 2a    */
#define XC_GGA_XC_SB98_2b  180 /* Schmider-Becke 98 parameterization 2b    */
#define XC_GGA_XC_SB98_2c  181 /* Schmider-Becke 98 parameterization 2c    */
#define XC_GGA_XC_HCTH_A    97 /* HCTH-A                                   */
#define XC_GGA_XC_B97_GGA1  96 /* Becke 97 GGA-1                           */
#define XC_GGA_XC_HCTH_P14  95 /* HCTH p=1/4                               */
#define XC_GGA_XC_HCTH_P76  94 /* HCTH p=7/6                               */
#define XC_GGA_XC_HCTH_407P 93 /* HCTH/407+                                */
#define XC_GGA_C_N12        80 /* N12 functional from Minnesota            */
#define XC_GGA_C_N12_SX     79 /* N12-SX functional from Minnesota         */
#define XC_GGA_XC_WB97     251 /* Chai and Head-Gordon                     */
#define XC_GGA_XC_WB97X    252 /* Chai and Head-Gordon                     */
#define XC_GGA_XC_WB97X_V  253 /* Mardirossian and Head-Gordon             */
#define XC_GGA_XC_WB97X_D  256 /* Chai and Head-Gordon                     */

static const FLOAT b97_params[][3][5] = {
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
  }, {   /* Becke 97 */
    {0.8094, 0.5073,  0.7481, 0.0, 0.0},                      /* X   */
    {0.1737, 2.3487, -2.4868, 0.0, 0.0},                      /* Css */
    {0.9454, 0.7471, -4.5961, 0.0, 0.0}                       /* Cab */
  }, {   /* Becke 97-1 */
    {0.789518, 0.573805,  0.660975, 0.0, 0.0},                /* X   */
    {0.0820011, 2.71681, -2.87103,  0.0, 0.0},                /* Css */
    {0.955689, 0.788552, -5.47869,  0.0, 0.0}                 /* Cab */
  }, {   /* Becke 97-2 */
    {0.827642,  0.0478400, 1.76125,  0.0, 0.0},               /* X   */
    {0.585808, -0.691682,  0.394796, 0.0, 0.0},               /* Css */
    {0.999849,  1.40626,  -7.44060,  0.0, 0.0}                /* Cab */
  }, {   /* Becke 97-D */
    {1.08662, -0.52127,  3.25429, 0.0, 0.0},                  /* X   */
    {0.22340, -1.56208,  1.94293, 0.0, 0.0},                  /* Css */
    {0.69041,  6.30270, -14.9712, 0.0, 0.0}                   /* Cab */
  }, {   /* Becke 97-K */
    {0.507863, 1.46873, -1.51301, 0.0, 0.0},                  /* X   */
    {0.12355,  2.65399, -3.20694, 0.0, 0.0},                  /* Css */
    {1.58613, -6.20977,  6.46106, 0.0, 0.0}                   /* Cab */
  }, {   /* Becke 97-3 */
    { 0.7334648,  0.2925270, 3.338789, -10.51158,  10.60907},  /* X   */
    { 0.5623649, -1.322980,  6.359191, -7.464002,   1.827082}, /* Css */
    { 1.133830,  -2.811967,  7.431302, -1.969342, -11.74423}   /* Cab */
  }, {   /* SB98-1a */
    { 0.845975,  0.228183,  0.749949, 0.0, 0.0},  /* X   */
    {-0.817637, -0.054676,  0.592163, 0.0, 0.0},  /* Css */
    { 0.975483,  0.398379, -3.73540,  0.0, 0.0}   /* Cab */
  }, {   /* SB98-1b */
    { 0.800103, -0.084192,  1.47742, 0.0, 0.0},  /* X   */
    { 1.44946,  -2.37073,   2.13564, 0.0, 0.0},  /* Css */
    { 0.977621,  0.931199, -4.76973, 0.0, 0.0}   /* Cab */
  }, {   /* SB98-1c */
    { 0.810936, 0.496090,  0.772385, 0.0, 0.0},  /* X   */
    { 0.262077, 2.12576,  -2.30465,  0.0, 0.0},  /* Css */
    { 0.939269, 0.898121, -4.91276,  0.0, 0.0}   /* Cab */
  }, {   /* SB98-2a */
    { 0.749200, 0.402322,  0.620779, 0.0, 0.0},  /* X   */
    { 1.26686,  1.67146,  -1.22565,  0.0, 0.0},  /* Css */
    { 0.964641, 0.050527, -3.01966,  0.0, 0.0}   /* Cab */
  }, {   /* SB98-2b */
    { 0.770587, 0.180767,  0.955246, 0.0, 0.0},  /* X   */
    { 0.170473, 1.24051,  -0.862711, 0.0, 0.0},  /* Css */
    { 0.965362, 0.863300, -4.61778,  0.0, 0.0}   /* Cab */
  }, {   /* SB98-2c */
    { 0.790194, 0.400271,  0.832857, 0.0, 0.0},  /* X   */
    {-0.120163, 2.82332,  -2.59412,  0.0, 0.0},  /* Css */
    { 0.934715, 1.14105,  -5.33398,  0.0, 0.0}   /* Cab */
  }, {   /* HCTH-A  */
    { 1.09878,  -2.51173,   0.0156233, 0.0,     0.0},  /* X   */
    { 0.0136823, 0.268920, -0.550769,  1.03947, 0.0},  /* Css */
    { 0.836897,  1.72051,  -2.78498,  -4.57504, 0.0}   /* Cab */
  }, {   /* B97 GGA-1  */
    { 1.1068, -0.8765,    4.2639, 0.0, 0.0},  /* X   */
    { 0.4883, -2.117,    2.3235,  0.0, 0.0},  /* Css */
    { 0.7961,  5.7060, -14.9820,  0.0, 0.0}   /* Cab */
  }, {   /* HCTH p=1/4  */
    { 1.03161,  -0.360781,   3.51994, -4.95944,  2.41165},  /* X   */
    { 2.82414,   0.0318843, -1.78512,  2.39795, -0.876909}, /* Css */
    { 0.0821827, 4.56466,  -13.5529,  13.3820,  -3.17493}   /* Cab */
  }, {   /* HCTH p=7/6  */
    { 1.16525,  -0.583033, 2.51769,   3.81278,   -5.45906}, /* X   */
    {-3.92143,  -1.10098, -0.0914050, -0.859723, 2.07184},  /* Css */
    { 0.192949, -5.73335, 50.8757,   135.475,  101.268}     /* Cab */
  }, {   /* HCTH 407p  */
    { 1.08018, -0.4117,   2.4368,   1.3890, -1.3529},  /* X   */
    { 0.80302, -1.0479,   4.9807, -12.890,   9.6446},  /* Css */
    { 0.73604,  3.0270, -10.075,   20.611, -29.418}    /* Cab */
  }, {   /* N12  */
    { 0.0,          0.0,          0.0,          0.0,          0.0},          /* X   */
    { 1.00000e+00, -5.53170e+00,  3.07958e+01, -5.64196e+01,  3.21250e+01},  /* Css; coefficients flipped in original paper! */
    { 1.00000e+00,  3.24511e+00, -2.52893e+01,  1.44407e+01,  1.96870e+01}   /* Cab; coefficients flipped in original paper! */
  }, {   /* N12-SX  */
    { 0.0,          0.0,          0.0,          0.0,          0.0},          /* X   */ 
    { 2.63373e+00, -1.05450e+00, -7.29853e-01,  4.94024e+00, -7.31760e+00},  /* Css; coefficients flipped in original paper! */
    { 8.33615e-01,  3.24128e+00, -1.06407e+01, -1.60471e+01,  2.51047e+01}   /* Cab; coefficients flipped in original paper! */
  }, {   /* wB97  */
    { 1.00000e+00,  1.13116e+00, -2.74915e+00,  1.20900e+01, -5.71642e+00},  /* X   */
    { 1.00000e+00, -2.55352e+00,  1.18926e+01, -2.69452e+01,  1.70927e+01},  /* Css */
    { 1.00000e+00,  3.99051e+00, -1.70066e+01,  1.07292e+00,  8.88211e+00}   /* Cab */
  }, {   /* wB97X  */
    { 8.42294e-01,  7.26479e-01,  1.04760e+00, -5.70635e+00,  1.32794e+01},  /* X   */
    { 1.00000e+00, -4.33879e+00,  1.82308e+01, -3.17430e+01,  1.72901e+01},  /* Css */
    { 1.00000e+00,  2.37031e+00, -1.13995e+01,  6.58405e+00, -3.78132e+00}   /* Cab */
  }, {   /* wB97X-V */
    { 0.833,        0.603,        1.194,        0.0,          0.0        },  /* X   */
    { 0.556,       -0.257,        0.0,          0.0,          0.0        },  /* Css */
    { 1.219,       -1.850,        0.0,          0.0,          0.0        }   /* Cab */
  }, {   /* wB97X-D */
    { 7.77964e-01,  6.61160e-01,  5.74541e-01, -5.25671e+00,  1.16386e+01},  /* X   */
    { 1.00000e+00, -6.90539e+00,  3.13343e+01, -5.10533e+01,  2.64423e+01},  /* Css */
    { 1.00000e+00,  1.79413e+00, -1.20477e+01,  1.40847e+01, -8.50809e+00}   /* Cab */
  }
};

typedef struct{
  const FLOAT (*cc)[5];
} gga_xc_b97_params;


static void 
gga_xc_b97_init(XC(func_type) *p)
{
  gga_xc_b97_params *params;

  assert(p != NULL);

  p->n_func_aux  = 2;
  p->func_aux    = (XC(func_type) **) malloc(2*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));
  p->func_aux[1] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_X,    XC_POLARIZED);
  XC(func_init)(p->func_aux[1], XC_LDA_C_PW, XC_POLARIZED);

  assert(p->params == NULL);
  p->params = malloc(sizeof(gga_xc_b97_params));
  params = (gga_xc_b97_params *)(p->params);

  switch(p->info->number){
  case XC_GGA_XC_HCTH_93:   p->func =  0;  break;
  case XC_GGA_XC_HCTH_120:  p->func =  1;  break;
  case XC_GGA_XC_HCTH_147:  p->func =  2;  break;
  case XC_GGA_XC_HCTH_407:  p->func =  3;  break;
  case XC_GGA_XC_B97:       p->func =  4;  break;
  case XC_GGA_XC_B97_1:     p->func =  5;  break;
  case XC_GGA_XC_B97_2:     p->func =  6;  break;
  case XC_GGA_XC_B97_D:     p->func =  7;  break;
  case XC_GGA_XC_B97_K:     p->func =  8;  break;
  case XC_GGA_XC_B97_3:     p->func =  9;  break;
  case XC_GGA_XC_SB98_1a:   p->func = 10;  break;
  case XC_GGA_XC_SB98_1b:   p->func = 11;  break;
  case XC_GGA_XC_SB98_1c:   p->func = 12;  break;
  case XC_GGA_XC_SB98_2a:   p->func = 13;  break;
  case XC_GGA_XC_SB98_2b:   p->func = 14;  break;
  case XC_GGA_XC_SB98_2c:   p->func = 15;  break;
  case XC_GGA_XC_HCTH_A:    p->func = 16;  break;
  case XC_GGA_XC_B97_GGA1:  p->func = 17;  break;
  case XC_GGA_XC_HCTH_P14:  p->func = 18;  break;
  case XC_GGA_XC_HCTH_P76:  p->func = 19;  break;
  case XC_GGA_XC_HCTH_407P: p->func = 20;  break;
  case XC_GGA_C_N12:        p->func = 21;  break;
  case XC_GGA_C_N12_SX:     p->func = 22;  break;
  case XC_GGA_XC_WB97:
    p->func = 23;
    XC(lda_x_set_params)(p->func_aux[0], 4.0/3.0, XC_NON_RELATIVISTIC, 0.4);
    break;
  case XC_GGA_XC_WB97X:
    p->func = 24;
    XC(lda_x_set_params)(p->func_aux[0], 4.0/3.0, XC_NON_RELATIVISTIC, 0.3);
    break;
  case XC_GGA_XC_WB97X_V:
    p->func = 25;
    XC(lda_x_set_params)(p->func_aux[0], 4.0/3.0, XC_NON_RELATIVISTIC, 0.3);
    break;
  case XC_GGA_XC_WB97X_D:
    p->func = 26;
    XC(lda_x_set_params)(p->func_aux[0], 4.0/3.0, XC_NON_RELATIVISTIC, 0.2);
    break;
  default:
    fprintf(stderr, "Internal error in gga_b97\n");
    exit(1);
    break;
  }

  params->cc = b97_params[p->func];
}


void 
XC(mgga_b97_func_g)(const FLOAT *cc, FLOAT gamma, FLOAT s, int order, FLOAT *g, FLOAT *dgds, FLOAT *d2gds2)
{
  FLOAT s2, dd, x, dxds, d2xds2, dgdx, d2gdx2;

  s2 = s*s;
  dd = 1.0 + gamma*s2;
  x  = gamma * s2/dd;

  *g = cc[0] + x*(cc[1] + x*(cc[2] + x*(cc[3] + x*cc[4])));

  if(order < 1) return;

  dxds  = gamma * 2.0*s/(dd*dd);
  dgdx  = cc[1] + x*(2.0*cc[2] + x*(3.0*cc[3] + x*4.0*cc[4]));
  *dgds = dgdx*dxds;

  if(order < 2) return;
  
  d2gdx2  = 2.0*cc[2] + x*(6.0*cc[3] + x*12.0*cc[4]);
  d2xds2  = 2.0*gamma*(1.0 - 3.0*gamma*s2)/(dd*dd*dd);
  *d2gds2 = d2gdx2*dxds*dxds + dgdx*d2xds2;
}


static inline void
func(const XC(func_type) *p, XC(gga_work_c_t) *r)
{
  static const FLOAT sign[2] = {1.0, -1.0};
  const FLOAT gamma[3] = {0.004, 0.2, 0.006}; /* Xs, Css, Cab */

  XC(lda_work_t) lda_pw[3], lda_x[3];
  const gga_xc_b97_params *params;
  FLOAT x_avg, aux, aux12;
  FLOAT fx, dfxdx, d2fxdx2, fcpar, dfcpardx, d2fcpardx2, fcper, dfcperdx, d2fcperdx2;
  FLOAT dx_avgdxs[2], d2x_avgdxs2[3];
  int is, js;
 
  params = (gga_xc_b97_params *)(p->params);

  /* first we get the parallel and perpendicular LDAs */
  XC(lda_stoll) (p->func_aux[1], XC(lda_c_pw_func), r->dens, r->zeta, r->order, lda_pw);
  XC(lda_stoll) (p->func_aux[0], XC(lda_x_func),    r->dens, r->zeta, r->order, lda_x);

  /* initialize to zero */
  r->f = 0.0;
  if(r->order >= 1){
    r->dfdrs = r->dfdz = r->dfdxs[0] = r->dfdxs[1] = r->dfdxt = 0.0;
  }
  if(r->order >= 2){
    r->d2fdrs2 = r->d2fdrsz = r->d2fdrsxt = r->d2fdrsxs[0] = r->d2fdrsxs[1] = 0.0;
    r->d2fdz2 = r->d2fdzxt = r->d2fdzxs[0] = r->d2fdzxs[1] = r->d2fdxt2 = 0.0;
    r->d2fdxtxs[0] = r->d2fdxtxs[1] = r->d2fdxs2[0] = r->d2fdxs2[1] = r->d2fdxs2[2] = 0.0;
  }

  /* now we calculate the g functions for exchange and parallel correlation */
  for(is = 0; is < 2; is++){
    if(r->ds[is] < p->info->min_dens) continue;

    XC(mgga_b97_func_g)(params->cc[0], gamma[0], r->xs[is], r->order, &fx, &dfxdx, &d2fxdx2);
    XC(mgga_b97_func_g)(params->cc[1], gamma[1], r->xs[is], r->order, &fcpar, &dfcpardx, &d2fcpardx2);

    r->f += lda_x[is].zk*fx + lda_pw[is].zk*fcpar;

    if(r->order < 1) continue;

    r->dfdrs     += lda_x[is].dedrs*fx  + lda_pw[is].dedrs*fcpar;
    r->dfdz      += lda_x[is].dedz *fx  + lda_pw[is].dedz *fcpar;
    r->dfdxs[is] += lda_x[is].zk*dfxdx  + lda_pw[is].zk*dfcpardx;

    if(r->order < 2) continue;
    
    js = (is == 0) ? 0 : 2;

    r->d2fdrs2      += lda_x[is].d2edrs2*fx  + lda_pw[is].d2edrs2*fcpar;
    r->d2fdrsz      += lda_x[is].d2edrsz*fx  + lda_pw[is].d2edrsz*fcpar;
    r->d2fdrsxs[is] += lda_x[is].dedrs*dfxdx + lda_pw[is].dedrs*dfcpardx;
    r->d2fdz2       += lda_x[is].d2edz2*fx   + lda_pw[is].d2edz2*fcpar;
    r->d2fdzxs[is]  += lda_x[is].dedz*dfxdx  + lda_pw[is].dedz*dfcpardx;
    r->d2fdxs2[js]  += lda_x[is].zk*d2fxdx2  + lda_pw[is].zk*d2fcpardx2;
  }

  /* and now we add the opposite-spin contribution */
  aux   = r->xs[0]*r->xs[0] + r->xs[1]*r->xs[1];
  aux12 = SQRT(aux);
  x_avg = aux12/M_SQRT2;

  XC(mgga_b97_func_g)(params->cc[2], gamma[2], x_avg, r->order, &fcper, &dfcperdx, &d2fcperdx2);

  r->f += lda_pw[2].zk*fcper;

  if(r->order < 1) return;

  dx_avgdxs[0] = r->xs[0]/(aux12*M_SQRT2);
  dx_avgdxs[1] = r->xs[1]/(aux12*M_SQRT2);

  r->dfdrs    += lda_pw[2].dedrs*fcper;
  r->dfdz     += lda_pw[2].dedz *fcper;
  r->dfdxs[0] += lda_pw[2].zk*dfcperdx*dx_avgdxs[0];
  r->dfdxs[1] += lda_pw[2].zk*dfcperdx*dx_avgdxs[1];

  if(r->order < 2) return;

  d2x_avgdxs2[0] =  r->xs[1]*r->xs[1]/(aux*aux12*M_SQRT2);
  d2x_avgdxs2[1] = -r->xs[0]*r->xs[1]/(aux*aux12*M_SQRT2);
  d2x_avgdxs2[2] =  r->xs[0]*r->xs[0]/(aux*aux12*M_SQRT2);

  r->d2fdrs2     += lda_pw[2].d2edrs2*fcper;
  r->d2fdrsz     += lda_pw[2].d2edrsz*fcper;
  r->d2fdrsxs[0] += lda_pw[2].dedrs*dfcperdx*dx_avgdxs[0];
  r->d2fdrsxs[1] += lda_pw[2].dedrs*dfcperdx*dx_avgdxs[1];
  r->d2fdz2      += lda_pw[2].d2edz2*fcper;
  r->d2fdzxs[0]  += lda_pw[2].dedz*dfcperdx*dx_avgdxs[0];
  r->d2fdzxs[1]  += lda_pw[2].dedz*dfcperdx*dx_avgdxs[1];
  r->d2fdxs2[0]  += lda_pw[2].zk*(d2fcperdx2*dx_avgdxs[0]*dx_avgdxs[0] + dfcperdx*d2x_avgdxs2[0]);
  r->d2fdxs2[1]  += lda_pw[2].zk*(d2fcperdx2*dx_avgdxs[0]*dx_avgdxs[1] + dfcperdx*d2x_avgdxs2[1]);
  r->d2fdxs2[2]  += lda_pw[2].zk*(d2fcperdx2*dx_avgdxs[1]*dx_avgdxs[1] + dfcperdx*d2x_avgdxs2[2]);
}


#include "work_gga_c.c"

const XC(func_info_type) XC(func_info_gga_xc_b97) = {
  XC_GGA_XC_B97,
  XC_EXCHANGE_CORRELATION,
  "Becke 97",
  XC_FAMILY_GGA,
  {&xc_ref_Becke1997_8554, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_1) = {
  XC_GGA_XC_B97_1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-1",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init,
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_2) = {
  XC_GGA_XC_B97_2,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-2",
  XC_FAMILY_GGA,
  {&xc_ref_Wilson2001_9233, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_d) = {
  XC_GGA_XC_B97_D,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-D",
  XC_FAMILY_GGA,
  {&xc_ref_Grimme2006_1787, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_k) = {
  XC_GGA_XC_B97_K,
  XC_EXCHANGE_CORRELATION,
  "Boese-Martin for Kinetics",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_3) = {
  XC_GGA_XC_B97_3,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-3",
  XC_FAMILY_GGA,
  {&xc_ref_Keal2005_121103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_93) = {
  XC_GGA_XC_HCTH_93,
  XC_EXCHANGE_CORRELATION,
  "HCTH/93",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_120) = {
  XC_GGA_XC_HCTH_120,
  XC_EXCHANGE_CORRELATION,
  "HCTH/120",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_147) = {
  XC_GGA_XC_HCTH_147,
  XC_EXCHANGE_CORRELATION,
  "HCTH/147",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_407) = {
  XC_GGA_XC_HCTH_407,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2001_5497, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_1a) = {
  XC_GGA_XC_SB98_1a,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1a)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_1b) = {
  XC_GGA_XC_SB98_1b,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1b)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_1c) = {
  XC_GGA_XC_SB98_1c,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1c)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_2a) = {
  XC_GGA_XC_SB98_2a,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2a)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_2b) = {
  XC_GGA_XC_SB98_2b,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2b)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_sb98_2c) = {
  XC_GGA_XC_SB98_2c,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2c)",
  XC_FAMILY_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_a) = {
  XC_GGA_XC_HCTH_A,
  XC_EXCHANGE_CORRELATION,
  "HCTH-A",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_b97_gga1) = {
  XC_GGA_XC_B97_GGA1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97 GGA-1",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2000_160, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_p14) = {
  XC_GGA_XC_HCTH_P14,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=1/4",
  XC_FAMILY_GGA, 
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_p76) = {
  XC_GGA_XC_HCTH_P76,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=7/6",
  XC_FAMILY_GGA,
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_hcth_407p) = {
  XC_GGA_XC_HCTH_407P,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407+",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2003_5965, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_c_n12) = {
  XC_GGA_C_N12,
  XC_CORRELATION,
  "Minnesota N12 functional",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2012_2310, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_xc_b97_init,
  NULL, NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_c_n12_sx) = {
  XC_GGA_C_N12_SX,
  XC_CORRELATION,
  "Minnesota N12-SX functional",
  XC_FAMILY_GGA,
  {&xc_ref_Peverati2012_16187, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  gga_xc_b97_init,
  NULL, NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_wb97) = {
  XC_GGA_XC_WB97,
  XC_EXCHANGE_CORRELATION,
  "wB97 worker function",
  XC_FAMILY_GGA,
  {&xc_ref_Chai2008_084106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_wb97x) = {
  XC_GGA_XC_WB97X,
  XC_EXCHANGE_CORRELATION,
  "wB97X worker function",
  XC_FAMILY_GGA,
  {&xc_ref_Chai2008_084106, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_wb97x_v) = {
  XC_GGA_XC_WB97X_V,
  XC_EXCHANGE_CORRELATION,
  "wB97X worker function",
  XC_FAMILY_GGA,
  {&xc_ref_Mardirossian2014_9904, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

const XC(func_info_type) XC(func_info_gga_xc_wb97x_d) = {
  XC_GGA_XC_WB97X_D,
  XC_EXCHANGE_CORRELATION,
  "wB97D worker function",
  XC_FAMILY_GGA,
  {&xc_ref_Chai2008_6615, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC,
  1e-23, 1e-32, 0.0, 1e-32,
  gga_xc_b97_init, 
  NULL,
  NULL,
  work_gga_c,
  NULL
};

