(* type: work_gga_x *)
(* prefix:
  gga_x_pbeint_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pbeint_params * )(p->params);
*)

mu := s -> params_a_muGE + (params_a_muPBE - params_a_muGE)* \
   params_a_alpha*s^2/(1.0 + params_a_alpha * s^2):

(* this is the gga_x_pbe expression *)
f0 := s -> 1.0 + params_a_kappa * (1.0 - params_a_kappa/(params_a_kappa + mu(s)*s^2)):
f  := x -> f0(X2S * x):
