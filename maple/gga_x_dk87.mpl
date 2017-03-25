(* type: work_gga_x *)
(* prefix:
  gga_x_dk87_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_dk87_params * )(p->params);
*)

betag := 7/(432*Pi*(6*Pi^2)^(1/3))/X_FACTOR_C:

f := x -> 1 + betag*x^2*(1 + params_a_a1*x^params_a_alpha)/(1 + params_a_b1*x^2):