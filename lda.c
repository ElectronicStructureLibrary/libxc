#include <stdlib.h>
#include <assert.h>

#include "util.h"

/* initialization */
void lda_init(lda_type *p, int functional, int nspin)
{
  /* sanity check */
  assert(functional == XC_LDA_C_WIGNER ||
	 functional == XC_LDA_C_RPA    ||
	 functional == XC_LDA_C_HL     ||
	 functional == XC_LDA_C_GL     ||
	 functional == XC_LDA_C_VWN    ||
	 functional == XC_LDA_C_PZ     ||
	 functional == XC_LDA_C_OB_PZ  ||
	 functional == XC_LDA_C_PW     ||
	 functional == XC_LDA_C_OB_PW  ||
	 functional == XC_LDA_C_LYP    ||
	 functional == XC_LDA_C_AMGB);
  
  assert(nspin==XC_UNPOLARIZED || nspin==XC_POLARIZED);
  p->nspin = nspin;
  
  /* initialize the functionals */
  switch(functional){
  case XC_LDA_C_WIGNER:
    lda_c_wigner_init(p);
    break;

  case XC_LDA_C_RPA:
    lda_c_rpa_init(p);
    break;

  case XC_LDA_C_HL:
    lda_c_hl_init(p);
    break;

  case XC_LDA_C_GL:
    lda_c_gl_init(p);
    break;

  case XC_LDA_C_VWN:
    lda_c_vwn_init(p);
    break;
    
  case XC_LDA_C_PZ:
    lda_c_pz_init(p);
    break;
    
  case XC_LDA_C_OB_PZ:
    lda_c_ob_pz_init(p);
    break;
    
  case XC_LDA_C_PW:
    lda_c_pw_init(p);
    break;
    
  case XC_LDA_C_OB_PW:
    lda_c_ob_pw_init(p);
    break;
    
  case XC_LDA_C_LYP:
    lda_c_lyp_init(p);
    break;

  case XC_LDA_C_AMGB:
    lda_c_amgb_init(p);
    break;
  }
}


void lda(lda_type *p, double *rho, double *ec, double *vc)
{
  double dens, zeta, rs;
  
  assert(p!=NULL);
  
  /* get the trace and the polarization of the density */
  if(p->func->number!=XC_LDA_X && p->func->number!=XC_LDA_C_XALPHA ){
    rho2dzeta(p->nspin, rho, &dens, &zeta);
    
    if(dens <= MIN_DENS){
      int i;
      *ec = 0.0;
      for(i=0; i<p->nspin; i++) vc[i] = 0.0;
      return;
    }
    
    rs = RS(dens); /* Wigner radius */
  }
  
  switch(p->func->number){
  case(XC_LDA_X):
    lda_x(p, rho, ec, vc);
    break;
    
  case XC_LDA_C_WIGNER: 
    lda_c_wigner(p, rs, ec, vc);
    break;
    
  case XC_LDA_C_RPA:
    lda_c_rpa(p, rs, ec, vc);
    break;
    
  case XC_LDA_C_HL:
  case XC_LDA_C_GL:
    lda_c_hl(p, rs, zeta, ec, vc);
    break;
    
  case XC_LDA_C_XALPHA:
    lda_c_xalpha(p, rho, ec, vc);
    break;
    
  case XC_LDA_C_VWN:
    lda_c_vwn(p, rs, zeta, ec, vc);
    break;
    
  case XC_LDA_C_PZ:
  case XC_LDA_C_OB_PZ:
    lda_c_pz(p, rs, zeta, ec, vc);
    break;
    
  case XC_LDA_C_PW:
  case XC_LDA_C_OB_PW:
    lda_c_pw(p, rs, zeta, ec, vc);
    break;
    
  case XC_LDA_C_LYP:
    lda_c_lyp(p, rs, ec, vc);
    break;
    
  case XC_LDA_C_AMGB:
    lda_c_amgb(p, rho, ec, vc);
    break;
  }
}
