(* type: work_gga_x *)
(* prefix:
  gga_x_b86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b86_params * )(p->params);
*)

f := x -> 1.0 + paramsbeta*x^2/(1.0 + paramsgamma*x^2)^paramsomega:
