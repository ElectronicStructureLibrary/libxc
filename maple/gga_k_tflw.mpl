(* type: work_gga_x *)
(* prefix:
  gga_k_tflw_params *params;
 
  assert(p->params != NULL);
  params = (gga_k_tflw_params * )(p->params);
*)

f := x -> params_a_gamma + (params_a_lambda/8.0)*x^2/K_FACTOR_C: