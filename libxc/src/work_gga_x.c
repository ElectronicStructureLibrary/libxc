/************************************************************************
  This file is to be included in GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3), this
  routine performs the necessary conversions between a functional of s
  and of rho.
************************************************************************/

#ifndef HEADER
#  define HEADER 1
#endif

static void 
work_gga_x(void *p_, double *rho, double *sigma,
	   double *e, double *vrho, double *vsigma)
{
  xc_gga_type *p = p_;

  double sfact, dens;
  int is;

  *e   = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact     = 1.0;
    vsigma[1] = 0.0; /* there are no cross terms in these functionals */
  }else
    sfact     = 2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    double gdm, ds, rho13;
    double x, f, dfdx, ldfdx;
    int js = is==0 ? 0 : 2;

    vrho[is]   = 0.0;
    vsigma[js] = 0.0;
    if(rho[is] < MIN_DENS) continue;

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;
  
    ds    = rho[is]/sfact;
    rho13 = pow(ds, 1.0/3.0);
    x     = gdm/(ds*rho13);

#if   HEADER == 1
    func(p, x, &f, &dfdx, &ldfdx);
#elif HEADER == 2
    /* this second header is useful for functionals that depend
       explicitly both on s and on sigma */
    func(p, x, gdm*gdm, &f, &dfdx, &ldfdx, &(vsigma[js]));
#endif

    (*e) += -sfact*X_FACTOR_C*(ds*rho13)*f;
      
    vrho[is] += -4.0/3.0*X_FACTOR_C*rho13*(f - dfdx*x);
    if(gdm>MIN_GRAD)
      vsigma[js] = -sfact*X_FACTOR_C*(ds*rho13)*(vsigma[js]/(sfact*sfact) + dfdx*x/(2.0*sigma[js]));
    else
      vsigma[js] = -X_FACTOR_C/(sfact*(ds*rho13))*ldfdx;
  }

  *e /= dens; /* we want energy per particle */
}
