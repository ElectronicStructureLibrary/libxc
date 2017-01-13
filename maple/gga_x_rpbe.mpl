(* type: work_gga_x *)
(* prefix:
  gga_x_rpbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_rpbe_params * )(p->params);
*)

f0 := s -> 1.0 + params_a_kappa * (1.0 - exp(-params_a_mu * s^2/params_a_kappa)):
f  := x -> f0(X2S * x):
