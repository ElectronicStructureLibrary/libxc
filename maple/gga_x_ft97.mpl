(* type: work_gga_c *)
(* prefix:
  gga_x_ft97_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_ft97_params * )(p->params);
*)

beta := (rs, z, xs) -> params_a_beta0 
  + params_a_beta1*sigma_spin(rs, z, xs)/(params_a_beta2 + sigma_spin(rs, z, xs)):

fx := (rs, z, xs) -> beta(rs, z, xs)*xs^2/sqrt(1 + 9*xs^2*beta(rs, z, xs)^2*arcsinh(xs^2)^2):

f0 := (rs, zeta, xt, xs0, xs1) -> 
  - (X_FACTOR_C + fx(rs,  zeta, xs0))*n_spin(rs,  zeta)^(4/3) 
  - (X_FACTOR_C + fx(rs, -zeta, xs1))*n_spin(rs, -zeta)^(4/3):

(* we want energy per particle *)
f := (rs, zeta, xt, xs0, xs1) -> f0(rs, zeta, xt, xs0, xs1)/n_total(rs):