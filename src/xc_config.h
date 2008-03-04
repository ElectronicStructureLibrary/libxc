#if SINGLE_PRECISION
#  define FLOAT float
#  define POW   powf
#  define ASINH asinhf
#  define XC(x) xc_s_ ## x
#  define XC_U(X) XC_S_ ## X
#else
#  define FLOAT double
#  define POW   pow
#  define ASINH asinh
#  define XC(x) xc_ ## x
#  define XC_U(X) XC_ ## X
#endif

#define XC_FC_FUNC2(a,b) FC_FUNC_(a,b) 
#define XC_FC_FUNC(a,b) XC_FC_FUNC2(XC(a), XC_U(b))
