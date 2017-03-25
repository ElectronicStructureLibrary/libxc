(* type: work_gga_x *)
(* prefix:
  gga_x_b88_params *params;

  assert(p->params != NULL);
  params = (gga_x_b88_params * )(p->params);
*)

# standard B88
$ifdef gga_x_b88_params
params_a_beta  := 0.0042:
params_a_gamma := 6:
$endif

f_b88 := x ->
  1 + (params_a_beta/X_FACTOR_C)*x^2/(1 + params_a_gamma*params_a_beta*x*arcsinh(x)):

f := x -> f_b88(x):
