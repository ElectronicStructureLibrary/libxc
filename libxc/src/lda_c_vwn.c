#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 LDA parametrization of Vosko, Wilk & Nusair
************************************************************************/

/* some constants         e_c^P      e_c^F      alpha_c */
typedef struct {
  double  A[3]; /* e_c^P, e_c^F, alpha_c */
  double  b[3];
  double  c[3];
  double x0[3];
  double  Q[3];
  double  fpp;
} vwn_consts_type;

/* These numbers are taken from the original reference, but divided by
     two to convert from Rydbergs to Hartrees */
static vwn_consts_type vwn_consts[2] = {
  /* VWN parametrization of the correlation energy */
  {
    { 0.0310907, 0.01554535,  0.0      }, /*  A */
    { 3.72744,   7.06042,     1.13107  }, /*  b */
    {12.9352,   18.0578,     13.0045   }, /*  c */
    {-0.10498,  -0.32500,    -0.0047584}, /* x0 */
    { 0.0,       0.0,         0.0      }, /*  Q */
    0.0 /* fpp */
  },
  /* VWN RPA */
  {
    { 0.0310907, 0.01554535,  0.0      }, /*  A */
    {13.0720,   20.1231,      1.06835  }, /*  b */
    {42.7198,  101.578,      11.4813   }, /*  c */
    {-0.409286, -0.743294,   -0.228344 }, /* x0 */
    { 0.0,       0.0,         0.0      }, /*  Q */
    0.0 /* fpp */
  }
};

/* initialization */
void init_vwn_constants(vwn_consts_type *X)
{
  int i;

  X->A[2] = -1.0/(6.0*M_PI*M_PI);
  for(i=0; i<3; i++){
    X->Q[i] = sqrt(4.0*X->c[i] - X->b[i]*X->b[i]);
  }
  X->fpp = 4.0/(9.0*(pow(2.0, 1.0/3.0) - 1));
}

static void lda_c_vwn_init(void *p_)
{
  init_vwn_constants(&vwn_consts[0]);
}

static void lda_c_vwn_rpa_init(void *p_)
{
  init_vwn_constants(&vwn_consts[1]);
}

/* Eq. (4.4) of [1] */
void ec_i(vwn_consts_type *X, int i, double x, double *ec, double *decdrs)
{
  double f1, f2, f3, fx, qx, xx0, tx, tt;
  
  f1  = 2.0*X->b[i]/X->Q[i];
  f2  = X->b[i]*X->x0[i]/(X->x0[i]*X->x0[i] + X->b[i]*X->x0[i] + X->c[i]);
  f3  = 2.0*(2.0*X->x0[i] + X->b[i])/X->Q[i];
  fx  = x*x + X->b[i]*x + X->c[i];  /* X(x) */
  qx  = atan(X->Q[i]/(2.0*x + X->b[i]));
  xx0 = x - X->x0[i];
  
  *ec = X->A[i]*(log(x*x/fx) + f1*qx - f2*(log(xx0*xx0/fx) + f3*qx));
  
  tx  = 2.0*x + X->b[i];
  tt  = tx*tx + X->Q[i]*X->Q[i];
  *decdrs = X->A[i]*(2.0/x - tx/fx - 4.0*X->b[i]/tt -
		     f2*(2.0/xx0 - tx/fx - 4.0*(2.0*X->x0[i] + X->b[i])/tt));
}

/* the functional */
void lda_c_vwn(void *p_, double *rho, double *ec, double *vc, double *fc)
{
  lda_type *p = (lda_type *)p_;

  double dens, zeta;
  double rs[2], ec_1, dec_1;
  int func;
  vwn_consts_type *X;

  assert(p!=NULL);

  rho2dzeta(p->nspin, rho, &dens, &zeta);

  func = p->func->number - XC_LDA_C_VWN;
  assert(func==0 || func==1);
  X = &vwn_consts[func];

  /* Wigner radius */
  rs[1] = RS(dens);     /* rs          */
  rs[0] = sqrt(rs[1]);  /* sqrt(rs)    */

  ec_i(X, 0, rs[0], &ec_1, &dec_1);
  
  if(p->nspin == XC_UNPOLARIZED){
    *ec   = ec_1;
    vc[0] = ec_1 - rs[0]/6.0*dec_1;
  }else{
    double ec_2, ec_3, dec_2, dec_3, fz, dfz, decdx, decdz;
    double t1, t2, z3, z4;
    
    ec_i(X, 1, rs[0], &ec_2, &dec_2);
    ec_i(X, 2, rs[0], &ec_3, &dec_3);
    
    fz  =  FZETA(zeta);
    dfz = DFZETA(zeta);
    
    z3 = pow(zeta, 3);
    z4 = z3*zeta;
    t1 = (fz/X->fpp)*(1.0 - z4);
    t2 = fz*z4;
    
    *ec   =  ec_1 +  ec_3*t1 + ( ec_2 -  ec_1)*t2;
    decdx = dec_1 + dec_3*t1 + (dec_2 - dec_1)*t2;
    decdz = (ec_3/X->fpp)*(dfz*(1.0 - z4) - 4.0*fz*z3) +
      (ec_2 - ec_1)*(dfz*z4 + 4.0*fz*z3);
    
    t1 = *ec - rs[0]/6.0*decdx;
    
    vc[0] = t1 + (1.0 - zeta)*decdz;
    vc[1] = t1 - (1.0 + zeta)*decdz;
  }
	
}

func_type func_lda_c_vwn = {
  XC_LDA_C_VWN,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair",
  "LDA",
  "S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)",
  lda_c_vwn_init,
  NULL,
  lda_c_vwn
};

func_type func_lda_c_vwn_rpa = {
  XC_LDA_C_VWN_RPA,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair (parametrization of the RPA energy)",
  "LDA",
  "S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)",
  lda_c_vwn_rpa_init,
  NULL,
  lda_c_vwn 
};

