#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Wigner's parametrization from the low density limit
************************************************************************/

static void lda_c_wigner(void *p_, double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a = -0.44, b = 7.8;
  double dens, zeta, rs;
  double t;
  
  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs    =  RS(dens); /* Wigner radius */
  t     =  b + rs;

  *ec   =  a/t;
  vc[0] = -a/(t*t);                         /* now contains d ec/d rs */
  
  vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
  if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to return something */
}

const xc_func_info_type func_info_lda_c_wigner = {
  XC_LDA_C_WIGNER,
  XC_CORRELATION,
  "Wigner",
  XC_FAMILY_LDA,
  "E.P. Wigner, Trans. Faraday Soc. 34, 678 (1938)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_wigner, /* lda  */
};


/************************************************************************
 Random Phase Approximation (RPA)
************************************************************************/

static void lda_c_rpa(void *p_, double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  static double a = 0.0311, b = -0.047, c = 0.009, d = -0.017;
  double dens, zeta, rs;
  double lrs;

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  rs    =  RS(dens); /* Wigner radius */
  lrs = log(rs);
  
  *ec   = a*lrs + b + c*rs*lrs + d*rs;
  vc[0] = a/rs + c*(lrs + 1.0) + d;         /* now contains d ec/d rs */
  
  vc[0] = *ec - rs/3.0*vc[0];               /* and now d ec/d rho */
  if(p->nspin==XC_POLARIZED) vc[1] = vc[0]; /* have to erturn something */
}

const xc_func_info_type func_info_lda_c_rpa = {
  XC_LDA_C_RPA,
  XC_CORRELATION,
  "Random Phase Approximation (RPA)",
  XC_FAMILY_LDA,
  "M. Gell-Mann and K.A. Brueckner, Phys. Rev. 106, 364 (1957)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_rpa,    /* lda  */
};


/************************************************************************
   L. Hedin and  B.I. Lundqvist
   O. Gunnarsson and B. I. Lundqvist
************************************************************************/

static void hl_f(int func, int i, double rs, double *ec, double *vc)
{
  static const 
    double r[2][2] = {{21.0,   21.0},     /* HL unpolarized only*/
		      {11.4,   15.9}};    /* GL */
  static const 
    double c[2][2] = {{ 0.0225, 0.0225},  /* HL unpolarized only */
		      { 0.0333, 0.0203}}; /* GL */
  
  double a, x, x2, x3;
  
  x   = rs/r[func][i];
  x2  = x*x;
  x3  = x2*x;
  
  a   = log(1.0 + 1.0/x);
  *ec = -c[func][i]*((1.0 + x3)*a - x2 + 0.5*x - 1.0/3.0);
  
  *vc = -c[func][i]*a;
}


static void lda_c_hl(void *p_, double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;

  double ecp, vcp;
  double dens, zeta, rs;
  int func = p->info->number - XC_LDA_C_HL;

  /* sanity check */
  assert(func==0 || func==1);
  
  rho2dzeta(p->nspin, rho, &dens, &zeta);
  rs = RS(dens); /* Wigner radius */  

  hl_f(func, 0, rs, &ecp, &vcp);
  
  if(p->nspin==XC_UNPOLARIZED){
    *ec   = ecp;
    vc[0] = vcp;
    
  }else{ /* XC_POLARIZED */
    double ecf, vcf, fz, dfz, t1;
    
    fz  =  FZETA(zeta);
    dfz = DFZETA(zeta);
    
    hl_f(func, 1, rs, &ecf, &vcf);
    
    *ec = ecp + (ecf - ecp)*fz;                  /* the energy    */
    t1  = vcp + (vcf - vcp)*fz;
    
    vc[0] = t1 + (ecf - ecp)*dfz*( 1.0 - zeta);  /* the potential */
    vc[1] = t1 + (ecf - ecp)*dfz*(-1.0 - zeta);
  }
}

const xc_func_info_type func_info_lda_c_hl = {
  XC_LDA_C_HL,
  XC_CORRELATION,
  "Hedin & Lundqvist",
  XC_FAMILY_LDA,
  /* can someone get me this paper, so I can find all coefficients? */
  "L. Hedin and B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_hl,     /* lda  */
};

const xc_func_info_type func_info_lda_c_gl = {
  XC_LDA_C_GL,
  XC_CORRELATION,
  "Gunnarson & Lundqvist",
  XC_FAMILY_LDA,
  "O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_hl,     /* lda  */
};


/************************************************************************
 Slater's Xalpha functional

    Exc = alpha Ex
************************************************************************/

/* This correlation functional, added to the exchange functional, produces
a total exchange and correlation functional, Exc, equal to 3/2 * alpha * Ex 
Setting alpha equal to one gives the *usual* Slater Xalpha functional,
whereas alpha equal to 2/3 just leaves the exchange functional unchanged */
static void lda_c_xalpha(void *p_, double *rho, double *ec, double *vc, double *fc)
{
  xc_lda_type *p = (xc_lda_type *)p_;
  double a = 1.5*p->alpha - 1.0;
  int i;

  lda_x(p, rho, ec, vc, fc);

  (*ec) *= a;

  for(i=0; i<p->nspin; i++) vc[i] *= a;

  if(fc != NULL)
    for(i=0; i<p->nspin*p->nspin; i++) fc[i] *= a;

}

const xc_func_info_type func_info_lda_c_xalpha = {
  XC_LDA_C_XALPHA,
  XC_CORRELATION,
  "Slater's Xalpha",
  XC_FAMILY_LDA,
  NULL,
  XC_PROVIDES_EXC | XC_PROVIDES_VXC | XC_PROVIDES_FXC,
  NULL,         /* init */
  NULL,         /* end  */
  lda_c_xalpha  /* lda */
};

void xc_lda_c_xalpha_init(xc_lda_type *p, int nspin, int dim, double alpha)
{
  p->alpha = alpha;
  xc_lda_x_init(p, nspin, dim, XC_NON_RELATIVISTIC);
  p->info = &func_info_lda_c_xalpha;
}
