(* type: work_mgga_x *)
(* prefix:
  mgga_x_ms_params *params;
 
  assert(pt->params != NULL);
  params = (mgga_x_ms_params * ) (pt->params);
*)

fa := a -> (1.0 - a^2)^3 / (1.0 + a^3 + params_a_b*a^6):
f0 := (p, c) -> 1.0 + params_a_kappa*(1.0 - params_a_kappa/(params_a_kappa + MU_GE*p + c)):

f := (rs, x, t, u) -> f0(X2S^2*x^2, 0.0) + \
  fa((t - x^2/8.0)/K_FACTOR_C)*(f0(X2S^2*x^2, params_a_c) - f0(X2S^2*x^2, 0.0)):
