(* type: work_gga_x *)
(* prefix:
  gga_x_b88_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_b88_params * )(p->params);
*)

$ifdef gga_x_b88_params
params_a_beta  := 0.0042:
params_a_gamma := 6.0:
$endif

b88_f := x -> 1.0 + params_a_beta/X_FACTOR_C*x^2/(1.0 + params_a_gamma*params_a_beta*x*arcsinh(x)):

f := x -> b88_f(x):
