(* type: work_gga_c *)
(* prefix:
  gga_x_kt_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_kt_params * )(p->params);
*)

fx := (rs, z, xs) -> 
   params_a_gamma*n_spin(rs, z)^(4.0/3.0)*xs^2/(n_spin(rs, z)^(4.0/3.0) + params_a_delta):

f0 := (rs, zeta, xt, xs0, xs1) -> 
  - (X_FACTOR_C - fx(rs,  zeta, xs0))*n_spin(rs,  zeta)^(4.0/3.0) 
  - (X_FACTOR_C - fx(rs, -zeta, xs1))*n_spin(rs, -zeta)^(4.0/3.0):

(* we want energy per particle *)
f := (rs, zeta, xt, xs0, xs1) -> f0(rs, zeta, xt, xs0, xs1)/n_total(rs):
