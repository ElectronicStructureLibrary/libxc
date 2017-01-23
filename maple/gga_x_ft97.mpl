(* type: work_gga_c *)
(* prefix:
  gga_x_ft97_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_ft97_params * )(p->params);
*)

beta := (rs, z, xs) -> params_a_beta0 
  + params_a_beta1*sigma_spin(rs, z, xs)/(params_a_beta2 + sigma_spin(rs, z, xs)):

fx := (rs, z, xs) -> beta(rs, z, xs)*xs^2/sqrt(1.0 + 9.0*xs^2*beta(rs, z, xs)^2*arcsinh(xs^2)^2):

f0 := (rs, zeta, xt, xs_0_, xs_1_) -> 
  - (X_FACTOR_C + fx(rs,  zeta, xs_0_))*n_spin(rs,  zeta)^(4.0/3.0) 
  - (X_FACTOR_C + fx(rs, -zeta, xs_1_))*n_spin(rs, -zeta)^(4.0/3.0):

(* we want energy per particle *)
f := (rs, zeta, xt, xs_0_, xs_1_) -> f0(rs, zeta, xt, xs_0_, xs_1_)/n_total(rs):