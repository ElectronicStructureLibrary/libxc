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

FLOAT XC(lambert_w)(FLOAT z)
{
  FLOAT w;
  int i;

  /* Sanity check - function is only defined for z >= -1/e */
  if(z + 1.0/M_E < -10*FLOAT_EPSILON) {
    fprintf(stderr,"Error - Lambert function called with argument z = %e.\n",z);
    exit(1);
  } else if(z < -1.0/M_E)
    /* Value of W(x) at x=-1/e is -1 */
    return -1.0;
  
  /* If z is small, go with the first terms of the power expansion
     (if z smaller than cube root of epsilon, z^4 will be zero to
     machine precision).
   */
  if(ABS(z) < CBRT(FLOAT_EPSILON))
    return z - z*z + 1.5*z*z*z;

  /* Initial guess. */
  if(z <= -0.3140862435046707) { /* Point where sqrt and Taylor polynomials match */
    /* Near the branching point: first terms in eqn (4.22) */
    w = SQRT(2.0*M_E*z + 2.0) - 1.0;
    
  } else if(z <= 1.149876485041417) { /* Point where Taylor and log expansion match */

    /* Taylor series around origin */
    w = z - z*z + 1.5*z*z*z;

  } else {
    /* Asymptotic expansion */
    FLOAT lnz = LOG(z);

    w = lnz - LOG(lnz);
  }

  /* Find result through iteration */
  for(i=0; i<10; i++){
    FLOAT expmw, dw;
    expmw = EXP(-w);
    
    /* Halley's equation, (5.9) in Corless et al */
    if( w != -1.0 )
      dw = - (w - z*expmw) / ( w + 1.0 - (w + 2.0)/(2.0*w + 2.0)*(w - z*expmw) );
    else
      dw = 0.0;

    w += dw;
    if(ABS(dw) < 10*FLOAT_EPSILON*(1.0 + ABS(w)))
      return w;
  }

  /* This should never happen! */
  fprintf(stderr, "%s\n%s\n", "lambert_w: iteration limit reached",
	  "Should never happen: execution aborted");
  exit(1);
}


struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

/* cheb_eval is defined in util.h */
