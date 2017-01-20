(* type: work_gga_x *)
(* prefix:
  gga_x_pw86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw86_params * )(p->params);
*)

f0 := s -> 1.0 + params_a_aa*s^2 + params_a_bb*s^4 + params_a_cc*s^6:
f  := x -> f0(X2S*x)^(1.0/15.0):