(* type: work_gga_x *)
(* prefix:
  gga_x_optx_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_optx_params * )(p->params);
*)

f := x-> params_a_a + params_a_b*(params_a_gamma*x^2/(1 + params_a_gamma*x^2))^2: