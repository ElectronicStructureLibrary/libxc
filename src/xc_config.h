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

#ifdef SINGLE_PRECISION
#  define FLOAT float
#  define POW   powf
#  define LOG   logf
#  define ABS   fabsf
#  define EXP   expf
#  define ERF   erff
#  define ERFC  erfcf
#  define SIN   sinf
#  define COS   cosf
#  define TAN   tanf
#  define ATAN  atanf
#  define ATAN2 atan2f
#  define ASINH asinhf

#ifdef HAVE_SQRTF
#  define SQRT  sqrtf
#else
#  define SQRT  sqrt
#endif

#ifdef HAVE_CBRTF
#define CBRT cbrtf
#elif defined(HAVE_CBRT)
#define CBRT cbrt
#else
#define CBRT(x) powf((x), 1.0/3.0)
#endif

#  define XC(x) xc_s_ ## x
#  define XC_U(X) XC_S_ ## X
#  define FLOAT_EPSILON FLT_EPSILON
#  define FLOAT_MIN FLT_MIN
#  define FLOAT_MAX FLT_MAX

#else
/* Double precision */

#  define FLOAT double
#  define POW   pow
#  define LOG   log
#  define ABS   fabs
#  define EXP   exp
#  define ERF   erf
#  define ERFC  erfc
#  define SIN   sin
#  define COS   cos
#  define TAN   tan
#  define ATAN  atan
#  define ATAN2 atan2
#  define ASINH asinh
#  define SQRT  sqrt

#ifdef HAVE_CBRT
#define CBRT cbrt
#else
#define CBRT(x) pow((x), 1.0/3.0)
#endif

#  define XC(x) xc_ ## x
#  define XC_U(X) XC_ ## X
#  define FLOAT_EPSILON DBL_EPSILON
#  define FLOAT_MIN DBL_MIN
#  define FLOAT_MAX DBL_MAX

#endif

#define XC_FC_FUNC2(a,b) FC_FUNC_(a,b)
#define XC_FC_FUNC(a,b) XC_FC_FUNC2(XC(a), XC_U(b))
