#include <stdio.h>
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

#if defined(LDA_SPEEDUP)
  /* This is for the new interpolation scheme */
  int i, k; int n = 600;
  double *x, *y, *y2, *y2d;
  double a, b, rpb, ea, z, vc[2], rho[2];

  x   = (double *)malloc(n*sizeof(double));
  y   = (double *)malloc(n*sizeof(double));
  y2  = (double *)malloc(n*sizeof(double));
  y2d = (double *)malloc(n*sizeof(double));

  a = 0.025;
  b = 0.0000001;

  x[0]  = -0.01; y[0]  = 0.0; y2[0] = 0.0; y2d[0] = 0.0;
  x[1]  = 0.0; y[1]  = 0.0; y2[1] = 0.0; y2d[1] = 0.0;

  (*p).acc = gsl_interp_accel_alloc ();

  if(p->nspin == XC_UNPOLARIZED){
    p->zeta_npoints = 1;

    (*p).energy = malloc(p->zeta_npoints*sizeof(gsl_spline));
    (*p).pot    = malloc(p->zeta_npoints*sizeof(gsl_spline));
    (*p).potd   = malloc(p->zeta_npoints*sizeof(gsl_spline));

    rpb = b; ea = exp(a);
    for (i = 2; i < n; i++) {
          x[i] = b* (exp(a*(i-1))-1.0);
          lda_work(p, &(x[i]), &(y[i]), &(y2[i]), NULL);}

    (*p).energy[0]     = gsl_spline_alloc (gsl_interp_linear, n);
    (*p).pot[0]        = gsl_spline_alloc (gsl_interp_linear, n);
    gsl_spline_init ((*p).energy[0], x, y, n);
    gsl_spline_init ((*p).pot[0], x, y2, n);

  }else{
    p->zeta_npoints = 11;

    (*p).energy = malloc(p->zeta_npoints*sizeof(gsl_spline));
    (*p).pot    = malloc(p->zeta_npoints*sizeof(gsl_spline));
    (*p).potd   = malloc(p->zeta_npoints*sizeof(gsl_spline));

    for(k = 0; k < p->zeta_npoints; k++){
      z = (1.0/(p->zeta_npoints-1))*k;
      rpb = b; ea = exp(a);
      for (i = 2; i < n; i++) {
            x[i] = b* (exp(a*(i-1))-1.0);
            rho[0] = 0.5*x[i]*(1.0+z);
            rho[1] = 0.5*x[i]*(1.0-z);
            lda_work(p, &(rho[0]), &(y[i]), &(vc[0]), NULL);
            y2[i] = vc[0]; y2d[i] = vc[1];
          }

      (*p).energy[k]     = gsl_spline_alloc (gsl_interp_linear, n);
      (*p).pot[k]        = gsl_spline_alloc (gsl_interp_linear, n);
      (*p).potd[k]       = gsl_spline_alloc (gsl_interp_linear, n);
      gsl_spline_init ((*p).energy[k], x, y, n);
      gsl_spline_init ((*p).pot[k], x, y2, n);
      gsl_spline_init ((*p).potd[k], x, y2d, n);
    }

  }
  free(x); free(y); free(y2); free(y2d);
#endif
}


void lda_work(lda_type *p, double *rho, double *ec, double *vc, double *fxc)
{
  double dens, zeta, rs;
  int got_fxc=0;

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
    lda_x(p, rho, ec, vc, fxc);
    got_fxc = 1;
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
    lda_c_xalpha(p, rho, ec, vc, fxc);
    got_fxc = 1;
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
    lda_c_pw(p, rs, dens, zeta, ec, vc, fxc);
    got_fxc = 1;
    break;
    
  case XC_LDA_C_LYP:
    lda_c_lyp(p, rs, ec, vc);
    break;
    
  case XC_LDA_C_AMGB:
    lda_c_amgb(p, rho, ec, vc);
    break;
  }
}

void lda(lda_type *p, double *rho, double *ec, double *vc)
{

#if defined(LDA_SPEEDUP)
  int i, j, k;
  double dens, zeta, rs;

  if(p->nspin == XC_UNPOLARIZED){
    *ec = gsl_spline_eval ( (*p).energy[0], rho[0],  (*p).acc );
    *vc = gsl_spline_eval ( (*p).pot[0], rho[0],  (*p).acc );
    return;}

  if(p->func->number == XC_LDA_X){
      *ec = gsl_spline_eval ( (*p).energy[0], rho[0],  (*p).acc);
       vc[0] = gsl_spline_eval ( (*p).pot[0], rho[0],  (*p).acc);
      *ec += gsl_spline_eval ( (*p).energy[0], rho[1], (*p).acc);
       vc[1] = gsl_spline_eval ( (*p).pot[0], rho[1],  (*p).acc);
       return;}

  rho2dzeta(p->nspin, rho, &dens, &zeta);
  if(zeta>0.0){ j=0; k=1; }else{ j=1; k=0; zeta=-zeta; };
  i = (int)zeta*(p->zeta_npoints -1);

  if( zeta == 0.0){
     *ec = gsl_spline_eval( (*p).energy[i], dens, (*p).acc);
      vc[0] = gsl_spline_eval ( (*p).pot[i], dens, (*p).acc);
      vc[1] = vc[0];
      return;}

  if( i == p->zeta_npoints-1){
     *ec = gsl_spline_eval( (*p).energy[i], dens, (*p).acc);
      vc[j] = gsl_spline_eval ( (*p).pot[i], dens, (*p).acc);
      vc[k] = gsl_spline_eval ( (*p).potd[i], dens, (*p).acc);
      return;}

  double ec1, ec2, vcu1, vcu2, vcd1, vcd2, dr;
  dr = zeta - i*(1.0/(p->zeta_npoints -1));
  ec1  = gsl_spline_eval( (*p).energy[i], dens, (*p).acc);
  vcu1 = gsl_spline_eval( (*p).pot[i],    dens, (*p).acc);
  vcd1 = gsl_spline_eval( (*p).potd[i],   dens, (*p).acc);
  ec2  = gsl_spline_eval( (*p).energy[i+1], dens, (*p).acc);
  vcu2 = gsl_spline_eval( (*p).pot[i+1],    dens, (*p).acc);
  vcd2 = gsl_spline_eval( (*p).potd[i+1],   dens, (*p).acc);
  *ec = ec1 + dr*(ec2-ec1)*(p->zeta_npoints -1);
  vc[j] = vcu1 + dr*(vcu2-vcu1)*(p->zeta_npoints -1);
  vc[k] = vcd1 + dr*(vcd2-vcd1)*(p->zeta_npoints -1);

#else
  lda_work(p, rho, ec, vc, NULL);
  return;
#endif

}

void lda_fxc(lda_type *p, double *rho, double *fxc)
{
  if(p->func->number == XC_LDA_X ||
     p->func->number == XC_LDA_C_XALPHA ||
     p->func->number == XC_LDA_C_PW ||
     p->func->number == XC_LDA_C_OB_PW){

    double ec, vc[2];
    lda_work(p, rho, &ec, vc, fxc);

  }else{ /* get fxc through a numerical derivative */
    int i, j;
    double delta_rho = 1e-5;
    
    for(i=0; i<p->nspin; i++){
      double rho2[2], e, vc1[2], vc2[2];

      j = (i+1) % 2;

      rho2[i] = rho[i] + delta_rho;
      rho2[j] = rho[j];
      lda_work(p, rho2, &e, vc1, NULL);

      rho2[i] = rho[i] - delta_rho;
      lda_work(p, rho2, &e, vc2, NULL);

      fxc __(i,i) = (vc1[i] - vc2[i])/(2.0*delta_rho);
      if(p->nspin == XC_POLARIZED)
	fxc __(i,j) = (vc1[j] - vc2[j])/(2.0*delta_rho);
    }
  }
}


void lda_kxc(lda_type *p, double *rho, double *kxc)
{
  /* Kxc, this is a third order tensor with respect to the densities */

  int i, j, k;
  const double delta_rho = 1e-4;

  for(i=0; i < p->nspin; i++){
    for(j=0; j < p->nspin; j++){
      for(k=0; k < p->nspin; k++){
    

	double rho2[2], e, vc1[2], vc2[2], vc3[2];
	int n, der_dir, func_dir;

	/* This is a bit tricky, we have to calculate a third mixed
	   partial derivative, one is calculated analitically (Vxc)
	   and the other two numerically. */
	   
	/* Here we reorder the indexes so that the numerical
	   derivative is taken with respect to the same variable */

	if(i!=j) {
	  if(j==k){
	    func_dir=i;
	    der_dir=j;
	  } else {
	    func_dir=j;
	    der_dir=i;
	  }
	} else {
	  func_dir=k;
	  der_dir=j;
	}
 

	for(n=0;n< p->nspin ; n++) rho2[n]=rho[n];

	lda_work(p, rho , &e, vc2, NULL);

	rho2[der_dir]+=delta_rho;
	lda_work(p, rho2, &e, vc1, NULL);
	
	rho2[der_dir]-=2.0*delta_rho;
	lda_work(p, rho2, &e, vc3, NULL);

	kxc ___(i,j,k)=(vc1[func_dir] - 2.0*vc2[func_dir] + vc3[func_dir])/(delta_rho*delta_rho);
	
      }
    }
  }

}
