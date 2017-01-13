(* type: work_gga_x *)
(* prefix:
  gga_x_b86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b86_params * )(p->params);
*)

f := x -> 1.0 + params_a_beta*x^2/(1.0 + params_a_gamma*x^2)^params_a_omega:
