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
#include "util.h"

/*
  Lambert W function. 
  adapted from the Fortran code of Rickard Armiento

  Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), 
         Adv. in Comp. Math. 5(4):329-359. 
*/

double lambert_w(double z)
{
  double result;
  int i;

  /* If z too low, go with the first term of the power expansion, z */
  if(z < 1.0e-20)
    return z;

  /* Inital guess */
  if(fabs(z + 1.0/M_E) > 1.45){
    /* Asymptotic expansion at 0 and Inf */
    result = log(z);
    result = result - log(result);
  }else{
    /* Series expansion about -1/e to first order */
    result = sqrt(2.0*M_E*z + 2.0) - 1.0;
  }

  /* Find result through iteration */
  for(i=0; i<10; i++){
    double p, t;

    p = exp(result);
    t = result*p - z;
    if( result != -1.0 )
      t = t/(p*(result + 1.0) - 0.5*(result + 2.0)*t/(result + 1.0));
    else
      t = 0.0;

    result = result - t;
    if(fabs(t) < (2.48e-14)*(1.0 + fabs(result)))
      return result;
  }

  /* This should never happen! */
  fprintf(stderr, "%s\n%s\n", "lambert_w: iteration limit reached",
	  "Should never happen: execution aborted");
  exit(1);
}


/*
  BESSI0: modified Bessel function I0(x).  See A&S 9.8.1, 9.8.2.
  from nemo_3.0.7/src/kernel/misc/besselfunc.c
*/

double bessi0(double x)
{
  double t, tt, ti, u;
  
  t = ABS(x) / 3.75;
  tt = t * t;
  if (tt < 1.0) {
    u = 1.0 +
      tt * (3.5156229 +
	    tt * (3.0899424 +
		  tt * (1.2067492 +
			tt * (0.2659732 +
			      tt * (0.0360768 +
				    tt * 0.0045813)))));
    return (u);
  } else {
    ti = 1.0 / t;
    u = 0.39894228 +
      ti * (0.01328592 +
	    ti * (0.00225319 +
		  ti * (-0.00157565 +
			ti * (0.00916281 +
			      ti * (-0.02057706 +
                                          ti * (0.02635537 +
                                                ti * (-0.01647633 +
                                                      ti * 0.00392377)))))));
    return (u * exp(ABS(x)) / sqrt(ABS(x)));
  }
}

/*
  BESSI1: modified Bessel function I1(x).  See A&S 9.8.3, 9.8.4.
  from nemo_3.0.7/src/kernel/misc/besselfunc.c
*/

double bessi1(double x)
{
  double t, tt, ti, u;
  
  t = x / 3.75;
  tt = t * t;
  if (tt < 1.0) {
    u = 0.5 +
      tt * (0.87890594 +
	    tt * (0.51498869 +
		  tt * (0.15084934 +
			tt * (0.02658733 +
			      tt * (0.00301532 +
				    tt * 0.00032411)))));
    return (u * x);
  } else {
    if (t < 0.0){
      fprintf(stderr, "bessi1(%g): negative argument", x);
      exit(1);
    }

    ti = 1.0 / t;
    u = 0.39894228 +
      ti * (-0.03988024 +
	    ti * (-0.00362018 +
		  ti * (0.00163801 +
			ti * (-0.01031555 +
			      ti * (0.02282967 +
				    ti * (-0.02895312 +
					  ti * (0.01787654 +
						ti * -0.00420059)))))));
    return (u * exp(ABS(x)) / sqrt(ABS(x)));
  }
}

/*
  BESSK0: modified Bessel function K0(x).  See A&S 9.8.5, 9.8.6.
  from nemo_3.0.7/src/kernel/misc/besselfunc.c
*/

double bessk0(double x)
{
  double t, tt, ti, u;
  
  if (x < 0.0){
    fprintf(stderr, "bessk0(%g): negative argument", x);
    exit(1);
  }

  t = x / 2.0;
  if (t < 1.0) {
    tt = t * t;
    u = -0.57721566 +
      tt * (0.42278420 +
	    tt * (0.23069756 +
		  tt * (0.03488590 +
			tt * (0.00262698 +
			      tt * (0.00010750 +
				    tt * 0.00000740)))));
    return (u - log(t) * bessi0(x));
  } else {
    ti = 1.0 / t;
    u = 1.25331414 +
      ti * (-0.07832358 +
	    ti * (0.02189568 +
		  ti * (-0.01062446 +
			ti * (0.00587872 +
			      ti * (-0.00251540 +
				    ti * 0.00053208)))));
    return (u * exp(-x) / sqrt(x));
  }
}


/*
  BESSK1: modified Bessel function K1(x).  See A&S 9.8.7, 9.8.8.
  from nemo_3.0.7/src/kernel/misc/besselfunc.c
*/

double bessk1(double x)
{
  double t, tt, ti, u;

  if (x < 0.0){
    fprintf(stderr, "bessk1(%g): negative argument", x);
    exit(1);
  }

  t = x / 2.0;
  if (t < 1.0) {
    tt = t * t;
    u = 1.0 +
      tt * (0.15443144 +
	    tt * (-0.67278579 +
		  tt * (-0.18156897 +
			tt * (-0.01919402 +
			      tt * (-0.00110404 +
				    tt * -0.00004686)))));
    return (u / x + log(t) * bessi1(x));
  } else {
    ti = 1.0 / t;
    u = 1.25331414 +
      ti * (0.23498619 +
	    ti * (-0.03655620 +
		  ti * (0.01504268 +
			ti * (-0.00780353 +
			      ti * (0.00325614 +
				    ti * -0.00068245)))));
    return (u * exp(-x) / sqrt(x));
  }
}

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

static inline double
cheb_eval(const double x, const double *cs, const int N)
{
  int i;
  double twox, b0, b1, b2;

  b1 = 0.0;
  b0 = 0.0;

  twox = 2.0*x;
  for(i=N-1; i>=0; i--){
    b2 = b1;
    b1 = b0;
    b0 = twox*b1 - b2 + cs[i];
  }

  return 0.5*(b0 - b2);
}


/* expint_E1 calculates the single precision exponential integral, E1(X), for
   positive single precision argument X and the Cauchy principal value
   for negative X.  If principal values are used everywhere, then, for
   all X,
   
     E1(X) = -Ei(-X)
   or
     Ei(X) = -E1(-X).

   Based on the SLATEC routine by W. Fullerton.
*/
static double AE11_data[39] = {
   0.121503239716065790,
  -0.065088778513550150,
   0.004897651357459670,
  -0.000649237843027216,
   0.000093840434587471,
   0.000000420236380882,
  -0.000008113374735904,
   0.000002804247688663,
   0.000000056487164441,
  -0.000000344809174450,
   0.000000058209273578,
   0.000000038711426349,
  -0.000000012453235014,
  -0.000000005118504888,
   0.000000002148771527,
   0.000000000868459898,
  -0.000000000343650105,
  -0.000000000179796603,
   0.000000000047442060,
   0.000000000040423282,
  -0.000000000003543928,
  -0.000000000008853444,
  -0.000000000000960151,
   0.000000000001692921,
   0.000000000000607990,
  -0.000000000000224338,
  -0.000000000000200327,
  -0.000000000000006246,
   0.000000000000045571,
   0.000000000000016383,
  -0.000000000000005561,
  -0.000000000000006074,
  -0.000000000000000862,
   0.000000000000001223,
   0.000000000000000716,
  -0.000000000000000024,
  -0.000000000000000201,
  -0.000000000000000082,
   0.000000000000000017
};

static double AE12_data[25] = {
   0.582417495134726740,
  -0.158348850905782750,
  -0.006764275590323141,
   0.005125843950185725,
   0.000435232492169391,
  -0.000143613366305483,
  -0.000041801320556301,
  -0.000002713395758640,
   0.000001151381913647,
   0.000000420650022012,
   0.000000066581901391,
   0.000000000662143777,
  -0.000000002844104870,
  -0.000000000940724197,
  -0.000000000177476602,
  -0.000000000015830222,
   0.000000000002905732,
   0.000000000001769356,
   0.000000000000492735,
   0.000000000000093709,
   0.000000000000010707,
  -0.000000000000000537,
  -0.000000000000000716,
  -0.000000000000000244,
  -0.000000000000000058
};

static double E11_data[19] = {
  -16.11346165557149402600,
    7.79407277874268027690,
   -1.95540581886314195070,
    0.37337293866277945612,
   -0.05692503191092901938,
    0.00721107776966009185,
   -0.00078104901449841593,
    0.00007388093356262168,
   -0.00000620286187580820,
    0.00000046816002303176,
   -0.00000003209288853329,
    0.00000000201519974874,
   -0.00000000011673686816,
    0.00000000000627627066,
   -0.00000000000031481541,
    0.00000000000001479904,
   -0.00000000000000065457,
    0.00000000000000002733,
   -0.00000000000000000108
};

static double E12_data[16] = {
  -0.03739021479220279500,
   0.04272398606220957700,
  -0.13031820798497005440,
   0.01441912402469889073,
  -0.00134617078051068022,
   0.00010731029253063780,
  -0.00000742999951611943,
   0.00000045377325690753,
  -0.00000002476417211390,
   0.00000000122076581374,
  -0.00000000005485141480,
   0.00000000000226362142,
  -0.00000000000008635897,
   0.00000000000000306291,
  -0.00000000000000010148,
   0.00000000000000000315
};

static double AE13_data[25] = {
  -0.605773246640603460,
  -0.112535243483660900,
   0.013432266247902779,
  -0.001926845187381145,
   0.000309118337720603,
  -0.000053564132129618,
   0.000009827812880247,
  -0.000001885368984916,
   0.000000374943193568,
  -0.000000076823455870,
   0.000000016143270567,
  -0.000000003466802211,
   0.000000000758754209,
  -0.000000000168864333,
   0.000000000038145706,
  -0.000000000008733026,
   0.000000000002023672,
  -0.000000000000474132,
   0.000000000000112211,
  -0.000000000000026804,
   0.000000000000006457,
  -0.000000000000001568,
   0.000000000000000383,
  -0.000000000000000094,
   0.000000000000000023
};

static double AE14_data[26] = {
  -0.18929180007530170,
  -0.08648117855259871,
   0.00722410154374659,
  -0.00080975594575573,
   0.00010999134432661,
  -0.00001717332998937,
   0.00000298562751447,
  -0.00000056596491457,
   0.00000011526808397,
  -0.00000002495030440,
   0.00000000569232420,
  -0.00000000135995766,
   0.00000000033846628,
  -0.00000000008737853,
   0.00000000002331588,
  -0.00000000000641148,
   0.00000000000181224,
  -0.00000000000052538,
   0.00000000000015592,
  -0.00000000000004729,
   0.00000000000001463,
  -0.00000000000000461,
   0.00000000000000148,
  -0.00000000000000048,
   0.00000000000000016,
  -0.00000000000000005
};

#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

double expint_e1(double x){
  const double xmaxt = -GSL_LOG_DBL_MIN;      /* XMAXT = -LOG (R1MACH(1)) */
  const double xmax  = xmaxt - log(xmaxt);    /* XMAX = XMAXT - LOG(XMAXT) */

  double e1;

  if(x <= -10.0)
    e1 = exp(-x)/x * (1.+ cheb_eval(20./x + 1., AE11_data, 39));
  else if(x <= -4.0)
    e1 = exp(-x)/x * (1. + cheb_eval((40./x + 7.)/3., AE12_data, 25));
  else if(x <= -1.0)
    e1 = -log(fabs(x)) + cheb_eval((2.*x + 5.)/3., E11_data, 19);
  else if(x == 0.0)
    exit(1); /* put some appropriate error condition */
  else if(x <= 1.0)
    e1 = (-log(fabs(x)) - 0.6875 + x) + cheb_eval(x, E12_data, 16);
  else if(x <= 4.0)
    e1 = exp(-x)/x * (1. + cheb_eval((8./x - 5.)/3., AE13_data, 25));
  else if(x <= xmax)
    e1 = exp(-x)/x * (1. + cheb_eval(8./x - 1., AE14_data, 26));
  else
    exit(1);

  return e1;
}
