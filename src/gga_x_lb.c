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


#include "util.h"

/* Note: Do not forget to add a correlation (LDA) functional to the
   LB94.
   
   Note 2: The 160 value is hardcoded in xc.h and libxc_master.F90 to
   define XC_GGA_XC_LB to keep backwards compatibility.

*/
#define XC_GGA_X_LB  160 /* van Leeuwen & Baerends */
#define XC_GGA_X_LBM 182 /* van Leeuwen & Baerends modified*/

typedef struct{
  int    modified; /* shall we use a modified version */
  double threshold; /* when to start using the analytic form */
  double ip;        /* ionization potential of the species */
  double qtot;      /* total charge in the region */

  double aa;     /* the parameters of LB94 */
  double gamm;

  double alpha;
  double beta;
} XC(gga_x_lb_params);

/************************************************************************
  Calculates van Leeuwen Baerends functional
************************************************************************/

static void
gga_lb_init(XC(func_type) *p)
{
  XC(gga_x_lb_params) *params;

  assert(p->params == NULL);

  p->n_func_aux  = 1;
  p->func_aux    = (XC(func_type) **) malloc(1*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_LDA_X, p->nspin);

  p->params = malloc(sizeof(XC(gga_x_lb_params)));

  params = (XC(gga_x_lb_params) *) (p->params);
  switch(p->info->number){
  case XC_GGA_X_LB:
    params->alpha = 1.0;
    params->beta  = 0.05;
    break;
  case XC_GGA_X_LBM:
    params->alpha = 1.19;
    params->beta  = 0.01;
    break;
  }
}


void 
XC(gga_lb_modified)(const XC(func_type) *func, int np, const double *rho, const double *sigma, double r, double *vrho)
{
  int ip, is, is2;
  double ds, gdm, x, sfact;

  XC(gga_x_lb_params) *params;

  assert(func != NULL);

  assert(func->params != NULL);
  params = (XC(gga_x_lb_params) *) (func->params);

  XC(lda_vxc)(func->func_aux[0], np, rho, vrho);

  sfact = (func->nspin == XC_POLARIZED) ? 1.0 : 2.0;

  for(ip=0; ip<np; ip++){
    for(is=0; is<func->nspin; is++){
      is2 = 2*is;

      vrho[is] *= params->alpha;

      gdm    = max(SQRT(sigma[is2])/sfact, MIN_GRAD);
      ds     = rho[is]/sfact;

      if(params->modified == 0 || 
	 (ds > params->threshold && gdm > params->threshold)){
	double f;
	
	if(ds <= func->dens_threshold) continue;
	
	x =  gdm/POW(ds, 4.0/3.0);
	
	if(x < 300.0) /* the actual functional */	   
	  f = -params->beta*x*x/(1.0 + 3.0*params->beta*x*ASINH(params->gamm*x));
	else          /* asymptotic expansion */
	  f = -x/(3.0*LOG(2.0*params->gamm*x));

	vrho[is] += f * CBRT(ds);
	
      }else if(r > 0.0){
	/* the asymptotic expansion of LB94 */
	x = r + (3.0/params->aa)*
	  LOG(2.0*params->gamm * params->aa * 1.0 / CBRT(params->qtot));
	
	/* x = x + POW(qtot*EXP(-aa*r), 1.0/3.0)/(beta*aa*aa); */
	
	vrho[is] -= 1.0/x;
      }
    }
    /* increment pointers */
    rho   += func->n_rho;
    sigma += func->n_sigma;
    
    if(vrho != NULL)
      vrho   += func->n_vrho;

  } /* ip loop */
}


static void 
gga_x_lb(const XC(func_type) *p, int np, const double *rho, const double *sigma,
	 double *zk, double *vrho, double *vsigma,
	 double *v2rho2, double *v2rhosigma, double *v2sigma2,
	 double *v3rho3, double *v3rho2sigma, double *v3rhosigma2, double *v3sigma3)
{
  XC(gga_lb_modified)(p, np, rho, sigma, 0.0, vrho);
}


static const func_params_type ext_params[] = {
  {  0, "Modified: 0 (no) | 1 (yes)"},
  {0.0, "Ionization potential (a.u.)"},
  {1e-32, "Threshold"},
  {0.0, "Total charge (necessary to fix the asymptotics"}
};


static void 
set_ext_params(XC(func_type) *p, const double *ext_params)
{
  XC(gga_x_lb_params) *params;
  double ff;

  assert(p!=NULL && p->params!=NULL);
  params = (XC(gga_x_lb_params) *) (p->params);

  ff = (ext_params == NULL) ? p->info->ext_params[0].value : ext_params[0];
  params->modified  = (int)round(ff);
  ff = (ext_params == NULL) ? p->info->ext_params[1].value : ext_params[1];
  params->threshold = ff;
  ff = (ext_params == NULL) ? p->info->ext_params[2].value : ext_params[2];
  params->ip        = ff;
  ff = (ext_params == NULL) ? p->info->ext_params[3].value : ext_params[3];
  params->qtot      = ff;

  if(params->modified){
    params->aa   = (params->ip > 0.0) ? 2.0*SQRT(2.0*params->ip) : 0.5;
    params->gamm = CBRT(params->qtot)/(2.0*params->aa);
  }else{
    params->aa   = 0.5;
    params->gamm = 1.0;
  }
}


const XC(func_info_type) XC(func_info_gga_x_lb) = {
  XC_GGA_X_LB,
  XC_EXCHANGE,
  "van Leeuwen & Baerends",
  XC_FAMILY_GGA,
  {&xc_ref_vanLeeuwen1994_2421, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  1e-32,
  4, ext_params, set_ext_params,
  gga_lb_init, NULL,
  NULL, gga_x_lb, NULL
};


const XC(func_info_type) XC(func_info_gga_x_lbm) = {
  XC_GGA_X_LBM,
  XC_EXCHANGE,
  "van Leeuwen & Baerends modified",
  XC_FAMILY_GGA,
  {&xc_ref_Schipper2000_1344, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_VXC,
  1e-32,
  4, ext_params, set_ext_params,
  gga_lb_init, NULL,
  NULL, gga_x_lb, NULL
};

