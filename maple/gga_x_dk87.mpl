(* type: work_gga_x *)
(* prefix:
  gga_x_dk87_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_dk87_params * )(p->params);
*)

betag := 7.0/(432.0*Pi*(6.0*Pi^2)^(1.0/3.0))/X_FACTOR_C:

f := x -> 1.0 + betag*x^2*(1.0 + params_a_a1*x^params_a_alpha)/(1.0 + params_a_b1*x^2):