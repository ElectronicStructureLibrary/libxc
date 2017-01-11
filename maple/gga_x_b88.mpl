(* type: work_gga_x *)
(* prefix:
  gga_x_b88_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b88_params * )(p->params);
*)

f0 := x -> 1.0 + paramsbeta/X_FACTOR_C*x^2/(1.0 + paramsgamma*paramsbeta*x*arcsinh(x)):

f := x -> f0(x):
