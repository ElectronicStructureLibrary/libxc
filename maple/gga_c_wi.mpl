(* type: work_gga_c *)
(* prefix:
  gga_c_wi_params *params;

  assert(p->params != NULL);
  params = (gga_c_wi_params * )(p->params);
*)

f_num := xt -> params_a_a + params_a_b*xt^2*exp(-params_a_k*xt^2):
f_den := (rs, xt) -> params_a_c + rs*(1.0 + params_a_d*(4.0*Pi/3.0)^(1.0/3.0)*xt^(7.0/2.0)):

f := (rs, zeta, xt, xs_1_, xs_2_) -> f_num(xt)/f_den(rs, xt):
