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

#ifndef _LDA_H
#define _LDA_H

/* These are generic header files that are needed basically everywhere */
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "xc.h"

/* xc_config.h needs to be included to use FLOAT and related macros*/
#include "xc_config.h"

/* we include the references also */
#include "references.h"

#ifndef M_E
# define M_E            2.7182818284590452354   /* e */
#endif
#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif
#ifndef M_SQRT2
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

#ifdef _MSC_VER
#define strcasecmp  _stricmp
#define strncasecmp _strnicmp

double asinh (double x);
float  asinhf(float  x);
double erf(double);
double erfc(double);
#endif


#define M_SQRTPI        1.772453850905516027298167483341145182798L
#define M_SQRT3         1.732050807568877293527446341505872366943L
#define M_CBRT2         1.259921049894873164767210607278228350570L
#define M_CBRT3         1.442249570307408382321638310780109588392L
#define M_CBRT4         1.587401051968199474751705639272308260391L
#define M_CBRT5         1.709975946676696989353108872543860109868L
#define M_CBRT6         1.817120592832139658891211756327260502428L
#define M_CBRT7         1.912931182772389101199116839548760282862L
#define M_CBRT9         2.080083823051904114530056824357885386338L

/* Very useful macros */
#ifndef min
#define min(x,y)  ((x<y) ? (x) : (y))
#endif
#ifndef max
#define max(x,y)  ((x<y) ? (y) : (x))
#endif

/* some useful constants */
#define LOG_DBL_MIN   (LOG(DBL_MIN))
#define LOG_DBL_MAX   (LOG(DBL_MAX))
#define SQRT_DBL_EPSILON   (SQRT(DBL_EPSILON))

/* special functions */
#define Heaviside(x) (((x) >= 0) ? 1.0 : 0.0)
double LambertW(double z);
FLOAT XC(dilogarithm)(const FLOAT x);

/* we define this function here, so it can be properly inlined by all compilers */
static inline FLOAT
XC(cheb_eval)(const FLOAT x, const FLOAT *cs, const int N)
{
  int i;
  FLOAT twox, b0, b1, b2;

  b2 = b1 = b0 = 0.0;

  twox = 2.0*x;
  for(i=N-1; i>=0; i--){
    b2 = b1;
    b1 = b0;
    b0 = twox*b1 - b2 + cs[i];
  }

  return 0.5*(b0 - b2);
}

FLOAT XC(bessel_I0_scaled)(const FLOAT x);
FLOAT XC(bessel_I0)(const FLOAT x);
FLOAT XC(bessel_K0_scaled)(const FLOAT x);
FLOAT XC(bessel_K0)(const FLOAT x);
FLOAT XC(bessel_K1_scaled)(const FLOAT x);
FLOAT XC(bessel_K1)(const FLOAT x);

FLOAT XC(expint_e1_impl)(const FLOAT x, const int scale);
static inline FLOAT expint_e1(const FLOAT x)         { return  XC(expint_e1_impl)( x, 0); }
static inline FLOAT expint_e1_scaled(const FLOAT x)  { return  XC(expint_e1_impl)( x, 1); }
static inline FLOAT expint_Ei(const FLOAT x)         { return -XC(expint_e1_impl)(-x, 0); }
#define Ei(x) expint_Ei(x)
static inline FLOAT expint_Ei_scaled(const FLOAT x)  { return -XC(expint_e1_impl)(-x, 1); }

/* integration */
typedef void integr_fn(FLOAT *x, int n, void *ex);
FLOAT XC(integrate)(integr_fn func, void *ex, FLOAT a, FLOAT b);
void XC(rdqagse)(integr_fn f, void *ex, FLOAT *a, FLOAT *b, 
	     FLOAT *epsabs, FLOAT *epsrel, int *limit, FLOAT *result,
	     FLOAT *abserr, int *neval, int *ier, FLOAT *alist__,
	     FLOAT *blist, FLOAT *rlist, FLOAT *elist, int *iord, int *last);
  
typedef struct XC(functional_key_t) {
  char name[256];
  int  number;
} XC(functional_key_t);


#define M_C 137.0359996287515 /* speed of light */

#define RS_FACTOR      0.6203504908994000166680068120477781673508     /* (3/(4*Pi))^1/3        */
#define X_FACTOR_C     0.9305257363491000250020102180716672510262     /* 3/8*cur(3/pi)*4^(2/3) */
#define X_FACTOR_2D_C  1.504505556127350098528211870828726895584      /* 8/(3*sqrt(pi))        */
#define K_FACTOR_C     4.557799872345597137288163759599305358515      /* 3/10*(6*pi^2)^(2/3)   */
#define MU_GE          0.1234567901234567901234567901234567901235     /* 10/81                 */
#define X2S            0.1282782438530421943003109254455883701296     /* 1/(2*(6*pi^2)^(1/3))  */
#define X2S_2D         0.1410473958869390717370198628901931464610     /* 1/(2*(4*pi)^(1/2))    */
#define FZETAFACTOR    0.5198420997897463295344212145564567011405     /* 2^(4/3) - 2           */

#define RS(x)          (RS_FACTOR/CBRT(x))
#define FZETA(x)       ((POW(1.0 + (x),  4.0/3.0) + POW(1.0 - (x),  4.0/3.0) - 2.0)/FZETAFACTOR)
#define DFZETA(x)      ((CBRT(1.0 + (x)) - CBRT(1.0 - (x)))*(4.0/3.0)/FZETAFACTOR)
#define D2FZETA(x)     ((4.0/9.0)/FZETAFACTOR)* \
  (ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -2.0/3.0) + pow(1.0 - (x), -2.0/3.0)))
#define D3FZETA(x)     (-(8.0/27.0)/FZETAFACTOR)* \
  (ABS(x)==1.0 ? (FLT_MAX) : (pow(1.0 + (x), -5.0/3.0) - pow(1.0 - (x), -5.0/3.0)))

#define MIN_DENS             5.0e-13
#define MIN_GRAD             5.0e-13
#define MIN_TAU              5.0e-13
#define MIN_ZETA             5.0e-13

/* The following inlines confuse the xlc compiler */
void XC(rho2dzeta)(int nspin, const FLOAT *rho, FLOAT *d, FLOAT *zeta);
void XC(fast_fzeta)(const FLOAT x, const int nspin, const int order, FLOAT * fz);
void XC(mix_init)(XC(func_type) *p, int n_funcs, const int *funcs_id, const FLOAT *mix_coef);

/* LDAs */
void XC(lda_init)(XC(func_type) *p);
void XC(lda_end) (XC(func_type) *p);

typedef struct XC(lda_work_t) {
  int   order; /* to which order should I return the derivatives */
  FLOAT rs, z;

  FLOAT f;                                   /* energy per unit particle */
  FLOAT dfdrs, dfdz;                         /*  first derivatives of e  */
  FLOAT d2fdrs2, d2fdrsz, d2fdz2;            /* second derivatives of e  */
  FLOAT d3fdrs3, d3fdrs2z, d3fdrsz2, d3fdz3; /*  third derivatives of e  */
} XC(lda_work_t);

void XC(lda_fxc_fd)(const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *fxc);
void XC(lda_kxc_fd)(const XC(func_type) *p, int np, const FLOAT *rho, FLOAT *kxc);

/* the different possibilities for screening the interaction */
#define XC_RSF_ERF      0
#define XC_RSF_ERF_GAU  1
#define XC_RSF_YUKAWA   2

typedef void XC(lda_func_type) (const XC(func_type) *p, XC(lda_work_t) *r);

void XC(lda_x_attenuation_function_erf)(int order, FLOAT aa, FLOAT *f, FLOAT *df, FLOAT *d2f, FLOAT *d3f);
void XC(lda_x_attenuation_function_erf_gau)(int order, FLOAT aa, FLOAT *f, FLOAT *df, FLOAT *d2f, FLOAT *d3f);
void XC(lda_x_attenuation_function_yukawa)(int order, FLOAT aa, FLOAT *f, FLOAT *df, FLOAT *d2f, FLOAT *d3f);
void XC(lda_x_attenuation_function)(int interaction, int order, FLOAT aa, FLOAT *f, FLOAT *df, FLOAT *d2f, FLOAT *d3f);

/* direct access to the internal functions */
void XC(lda_x_func)     (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_x_erf_func) (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_hl_func)  (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_vwn_func) (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_pw_func)  (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_pz_func)  (const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_rc04_func)(const XC(func_type) *p, XC(lda_work_t) *r);
void XC(lda_c_2d_amgb_func)(const XC(func_type) *p, XC(lda_work_t) *r);

/* GGAs */
typedef struct XC(gga_work_x_t) {
  int   order; /* to which order should I return the derivatives */
  FLOAT x;

  FLOAT f;          /* enhancement factor       */
  FLOAT dfdx;       /* first derivatives of f  */
  FLOAT d2fdx2;     /* second derivatives of zk */
  FLOAT d3fdx3;
} XC(gga_work_x_t);

void work_gga_becke_init(XC(func_type) *p);

/* exchange enhancement factors: if you add one, please add it also to the util.c */
typedef void(*xc_gga_enhancement_t)(const XC(func_type) *, XC(gga_work_x_t) *r);
xc_gga_enhancement_t XC(get_gga_enhancement_factor)(int func_id);

void XC(gga_x_wc_enhance)   (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_pbe_enhance)  (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_pw91_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_rpbe_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_htbs_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_b86_enhance)  (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_b88_enhance)  (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_g96_enhance)  (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_pw86_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_airy_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_ak13_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_bayesian_enhance)(const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_bpccac_enhance)(const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_c09x_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_am05_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_dk87_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_herman_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_lg93_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_lv_rpw86_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_mpbe_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_optx_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_sogga11_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_ssb_sw_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);
void XC(gga_x_vmt_enhance) (const XC(func_type) *p, XC(gga_work_x_t) *r);

/* these functions are used in more than one functional */
void XC(lda_c_pw_g)(int func, int order, int k, FLOAT *rs, FLOAT *f, FLOAT *dfdrs, FLOAT *d2fdrs2, FLOAT *d3fdrs3);
void XC(beta_Hu_Langreth) (FLOAT r, int order, FLOAT *b, FLOAT *dbdr, FLOAT *d2bdr2);

typedef struct XC(gga_work_c_t) {
  int   order; /* to which order should I return the derivatives */

  FLOAT dens, ds[2], sigmat, sigmas[3];
  FLOAT rs, z, xt, xs[2];

  FLOAT f;

  FLOAT dfdrs, dfdz, dfdxt, dfdxs[2];
  FLOAT d2fdrs2, d2fdrsz, d2fdrsxt, d2fdrsxs[2], d2fdz2, 
    d2fdzxt, d2fdzxs[2], d2fdxt2, d2fdxtxs[2], d2fdxs2[3];

  FLOAT d3fdrs3, d3fdz3, d3fdxt3, d3fdxs3[4]; /* uuu, uud, udd, ddd */
  FLOAT d3fdrs2z, d3fdrs2xt, d3fdrs2xs[2];
  FLOAT d3fdrsz2, d3fdz2xt, d3fdz2xs[2];
  FLOAT d3fdrsxt2, d3fdzxt2, d3fdxt2xs[2];
  FLOAT d3fdrsxs2[3], d3fdzxs2[3],d3fdxtxs2[3];
  FLOAT d3fdrszxt, d3fdrszxs[2], d3fdrsxtxs[2], d3fdzxtxs[2];
} XC(gga_work_c_t);

void XC(gga_c_pw91_func)(const XC(func_type) *p, XC(gga_work_c_t) *r);
void XC(gga_c_pbe_func) (const XC(func_type) *p, XC(gga_work_c_t) *r);
void XC(gga_c_pbeloc_func) (const XC(func_type) *p, XC(gga_work_c_t) *r);
void XC(gga_c_regtpss_func) (const XC(func_type) *p, XC(gga_work_c_t) *r);
void XC(gga_c_scan_e0_func) (const XC(func_type) *p, XC(gga_work_c_t) *r);
void XC(gga_c_q2d_func) (const XC(func_type) *p, XC(gga_work_c_t) *r);

/* meta GGAs */
typedef struct XC(mgga_work_x_t) {
  int   order; /* to which order should I return the derivatives */
  FLOAT rs, zeta, x, t, u;

  FLOAT f;                                   /* enhancement factor       */
  FLOAT dfdrs, dfdx, dfdt, dfdu;             /* first derivatives of f  */
  FLOAT d2fdrs2, d2fdx2, d2fdt2, d2fdu2;     /* second derivatives of zk */
  FLOAT d2fdrsx, d2fdrst, d2fdrsu, d2fdxt, d2fdxu, d2fdtu;
} XC(mgga_work_x_t);

typedef struct XC(mgga_work_c_t) {
  int   order; /* to which order should I return the derivatives */

  FLOAT dens, ds[2], sigmat, sigmas[3];
  FLOAT rs, z, xt, xs[2], ts[2], us[2];

  FLOAT f;
  FLOAT dfdrs, dfdz, dfdxt, dfdxs[2], dfdts[2], dfdus[2];
  FLOAT d2fdrs2, d2fdrsz, d2fdrsxt, d2fdrsxs[2], d2fdrsts[2], d2fdrsus[2];
  FLOAT d2fdz2, d2fdzxt, d2fdzxs[2], d2fdzts[2], d2fdzus[2];
  FLOAT d2fdxt2, d2fdxtxs[2], d2fdxtts[2], d2fdxtus[2];
  FLOAT d2fdxs2[3], d2fdxsts[4], d2fdxsus[4];
  FLOAT d2fdts2[3], d2fdtsus[4];
  FLOAT d2fdus2[3];
  FLOAT d3fdrs3, d3fdrs2z, d3fdrsz2, d3fdrszxt, d3fdrszxs[2], d3fdrszts[2], d3fdrszus[2];
  FLOAT d3fdrs2xt, d3fdrsxt2, d3fdrsxtxs[2], d3fdrsxtts[2], d3fdrsxtus[2], d3fdrs2xs[2];
  FLOAT d3fdrsxs2[3], d3fdrsxsts[4], d3fdrsxsus[4], d3fdrs2ts[2], d3fdrsts2[3];
  FLOAT d3fdrstsus[4], d3fdrs2us[2], d3fdrsus2[3];
  FLOAT d3fdz3, d3fdz2xt, d3fdzxt2, d3fdzxtxs[2], d3fdzxtts[2], d3fdzxtus[2];
  FLOAT d3fdz2xs[2], d3fdzxs2[3], d3fdzxsts[4], d3fdzxsus[4], d3fdz2ts[2], d3fdzts2[3];
  FLOAT d3fdztsus[4], d3fdz2us[2], d3fdzus2[3];
  FLOAT d3fdxt3, d3fdxt2xs[2], d3fdxtxs2[3], d3fdxtxsts[4], d3fdxtxsus[4], d3fdxt2ts[2];
  FLOAT d3fdxtts2[3], d3fdxttsus[4], d3fdxt2us[2], d3fdxtus2[3];
  FLOAT d3fdxs3[4], d3fdxs2ts[6], d3fdxs2us[6], d3fdxsts2[6], d3fdxstsus[8], d3fdxsus2[6];
  FLOAT d3fdts3[4], d3fdts2us[6], d3fdtsus2[6], d3fdus3[4];
} XC(mgga_work_c_t);


void XC(mgga_x_scan_falpha)(int order, FLOAT a, FLOAT c1, FLOAT c2, FLOAT dd, FLOAT *f, FLOAT *dfda);
	 

/* now the routines to set the _internal_ parameters of several functionals */
void XC(gga_x_pw91_set_params)(XC(func_type) *p, FLOAT a, FLOAT b, FLOAT c, FLOAT d, FLOAT f, FLOAT alpha, FLOAT expo);
void XC(gga_x_pw91_set_params2)(XC(func_type) *p, FLOAT bt, FLOAT alpha, FLOAT expo);
void XC(gga_x_ssb_sw_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT C, FLOAT D, FLOAT E);
void XC(gga_x_hjs_set_params)(XC(func_type) *p, FLOAT omega);
void XC(gga_x_hjs_b88_v2_set_params)(XC(func_type) *p, FLOAT omega);
void XC(gga_x_b88_set_params)(XC(func_type) *p, FLOAT beta, FLOAT gamma);
void XC(gga_x_optx_set_params)(XC(func_type) *p, FLOAT a, FLOAT b, FLOAT gamma);
void XC(gga_x_wpbeh_set_params)(XC(func_type) *p, FLOAT omega);
void XC(gga_x_pbe_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT mu);
void XC(gga_x_pbeint_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT alpha, FLOAT muPBE, FLOAT muGE);
void XC(gga_x_ityh_set_params)(XC(func_type) *p, int func_id, FLOAT omega);
void XC(gga_x_b86_set_params)(XC(func_type) *p, FLOAT beta, FLOAT gamma, FLOAT omega);
void XC(gga_x_rpbe_set_params)(XC(func_type) *p, FLOAT kappa, FLOAT mu);
void XC(gga_x_sfat_set_params)(XC(func_type) *p, int func_id, FLOAT omega);
void XC(gga_x_kt_set_params)(XC(func_type) *p, FLOAT gamma, FLOAT delta);
void XC(gga_c_lyp_set_params)(XC(func_type) *p, FLOAT A, FLOAT B, FLOAT c, FLOAT d);
void XC(gga_c_pbe_set_params)(XC(func_type) *p, FLOAT beta);

void XC(mgga_x_tpss_set_params)(XC(func_type) *p, FLOAT b, FLOAT c, FLOAT e, FLOAT kappa, FLOAT mu, FLOAT BLOC_a, FLOAT BLOC_b);
void XC(mgga_c_tpss_set_params)(XC(func_type) *p, FLOAT beta, FLOAT d, FLOAT C0_0, FLOAT C0_1, FLOAT C0_2, FLOAT C0_3);
void XC(mgga_c_pkzb_set_params)(XC(func_type) *p, FLOAT beta, FLOAT d, FLOAT C0_0, FLOAT C0_1, FLOAT C0_2, FLOAT C0_3);
void XC(mgga_c_bc95_set_params)(XC(func_type) *p, FLOAT css, FLOAT copp);

/* useful MACROS */
#define DFRACTION(num, dnum, den, dden) \
  (((dnum)*(den) - (num)*(dden))/((den)*(den)))
#define D2FRACTION(num, dnum, d2num, den, dden, d2den) \
  ((2.0*(num)*(dden)*(dden) - 2.0*(den)*(dden)*(dnum) - (den)*(num)*(d2den) + (den)*(den)*(d2num))/((den)*(den)*(den)))
#define D3FRACTION(num, dnum, d2num, d3num, den, dden, d2den, d3den)	\
  ((-(num)*(6.0*(dden)*(dden)*(dden) - 6.0*(den)*(dden)*(d2den) + (den)*(den)*(d3den)) + \
    (den)*(6.0*(dden)*(dden)*(dnum) - 3.0*(den)*(dden)*(d2num) + (den)*(-3.0*(dnum)*(d2den) + (den)*(d3num))))/((den)*(den)*(den)*(den)))


#endif
