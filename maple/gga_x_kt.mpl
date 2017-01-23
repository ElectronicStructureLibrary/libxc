(* type: work_gga_c *)
(* prefix:
  gga_x_kt_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_kt_params * )(p->params);
*)

n  := rs -> (RS_FACTOR/rs)^3:
nu := (rs, z) -> (1.0 + z)*n(rs)/2.0:
nd := (rs, z) -> (1.0 - z)*n(rs)/2.0:

f0 := (rs, zeta, xt, xs_0_, xs_1_) -> 
  - X_FACTOR_C*(nu(rs, zeta)^(4.0/3.0) + nd(rs, zeta)^(4.0/3.0))
  + params_a_gamma*nu(rs, zeta)^(8.0/3.0)*xs_0_^2/(nu(rs, zeta)^(4.0/3.0) + params_a_delta)
  + params_a_gamma*nd(rs, zeta)^(8.0/3.0)*xs_1_^2/(nd(rs, zeta)^(4.0/3.0) + params_a_delta):

(* we want energy per particle *)
f := (rs, zeta, xt, xs_0_, xs_1_) -> f0(rs, zeta, xt, xs_0_, xs_1_)/n(rs):
