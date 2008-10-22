/*
 Copyright (C) 2006-2008 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/************************************************************************
  This file is to be included in meta GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3) and tau, this
  routine performs the necessary conversions between a functional of s and tau
  and of rho.
************************************************************************/

static void
work_mgga_c_init(void *p_)
{
  XC(mgga_type) *p = (XC(mgga_type) *)p_;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_init)(p->lda_aux, XC_LDA_C_PW, XC_POLARIZED);
}


static void
work_mgga_c_end(void *p_)
{
  XC(mgga_type) *p = (XC(mgga_type) *)p_;

  /* XC(lda_end)(p->lda_aux); */
  free(p->lda_aux);
}


static void 
work_mgga_c(const void *p_, const FLOAT *rho, const FLOAT *sigma, const FLOAT *tau,
	    FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vtau,
	    FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2, FLOAT *v2rhotau, FLOAT *v2tausigma, FLOAT *v2tau2)
{
  const XC(mgga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  FLOAT ds[2], sigmas[2], x[2], t[2], f_LDA[2], vrho_LDA[2];
  int is, order;

  order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  if(p->nspin == XC_UNPOLARIZED)
    ds[1] = rho[0]/2.0;

  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, rho13;
    FLOAT f, ltau, dfdx, dfdt, d2fdx2, d2fdxt, d2fdt2;
    int js = (is == 0) ? 0 : 2;

    dens  += rho[is];
    ds[is] = rho[is]/sfact;

    if(rho[is] < MIN_DENS) continue;

    sigmas[is] = max(MIN_GRAD*MIN_GRAD, sigma[js]/sfact2);
    gdm        = sqrt(sigmas[is]);
  
    rho13  = POW(ds[is], 1.0/3.0);
    x [is] = gdm/(ds[is]*rho13);
    
    ltau   = max(tau[is]/sfact, MIN_TAU);
    t [is] = ltau/(ds[is]*rho13*rho13);  /* tau/rho^(5/3) */

    dfdx  = d2fdx2 = 0.0;

    dfdt = 0.0;
    func_c_parallel(p, x[is], t[is], order, &f, &dfdx, &dfdt, &d2fdx2, &d2fdxt, &d2fdt2);

    { /* get parallel spin LDA energy */
      FLOAT tmp_rho[2], tmp_vrho[2];

      tmp_rho[0] = ds[is];
      tmp_rho[1] = 0.0;

      switch (order){
      case 0:
	XC(lda_exc)(p->lda_aux, tmp_rho, &(f_LDA[is]));
	break;
      case 1:
	XC(lda_vxc)(p->lda_aux, tmp_rho, &(f_LDA[is]), tmp_vrho);
	vrho_LDA[is] = tmp_vrho[0];
	break;
      case 2: /* to be implemented */
	break; 
      }
    }

    if(zk != NULL)
      *zk += sfact*ds[is]*f_LDA[is]*f;
 
    if(vrho != NULL){
      vrho[is]    = vrho_LDA[is]*f - f_LDA[is]*(4.0*dfdx*x[is] + 5.0*dfdt*t[is])/3.0;
      vtau[is]    = f_LDA[is]*dfdt/(rho13*rho13);

      vsigma[js]  = ds[is]*f_LDA[is]*dfdx*x[is]/(2.0*sfact*sigmas[is]);
    }

    if(v2rho2 != NULL){
      /* Missing terms here */
      exit(1);
    }
  }
  /* *zk /= dens; return; */  /* DEBUG */

  /* We are now missing the opposite-spin part */
  {
    FLOAT f_LDA_opp, vrho_LDA_opp[2];
    FLOAT f, dfdx, dfdt, d2fdx2, d2fdxt, d2fdt2;
    FLOAT xt, tt;

    switch (order){
    case 0:
      XC(lda_exc)(p->lda_aux, ds, &f_LDA_opp);
      break;
    case 1:
      XC(lda_vxc)(p->lda_aux, ds, &f_LDA_opp, vrho_LDA_opp);
      break;
    case 2: /* to be implemented */
      break; 
    }

    if(p->nspin == XC_POLARIZED){
      xt = tt = 0.0;
      for(is=0; is<p->nspin; is++)
	if(rho[is] > MIN_DENS){
	  xt += x[is]*x[is];
	  tt += t[is];
	}
      xt = sqrt(xt);
    }else{
      xt = sqrt(2.0)*x[0];
      tt =      2.0 *t[0];
    }

    func_c_opposite(p, xt, tt, order, &f, &dfdx, &dfdt, &d2fdx2, &d2fdxt, &d2fdt2);

    if(zk != NULL)
      *zk += dens*f_LDA_opp*f;
 
    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	int js = (is == 0) ? 0 : 2;

	if(rho[is] < MIN_DENS) continue;

	vrho[is]   += vrho_LDA_opp[is]*f - dens*f_LDA_opp*(4.0*dfdx*x[is]*x[is]/xt + 5.0*dfdt*t[is])/(3.0*ds[is]);
	vtau[is]   += f_LDA_opp*dfdt*dens/POW(ds[is], 5.0/3.0);
	vsigma[js] += dens*f_LDA_opp*dfdx*x[is]*x[is]/(2.0*xt*sfact*sigmas[is]);
      }
    }
  }

  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}
