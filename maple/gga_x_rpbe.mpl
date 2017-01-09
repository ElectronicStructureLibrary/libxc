(* type: work_gga_x *)
(* prefix:
  gga_x_rpbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_rpbe_params * )(p->params);
*)

f0 := s -> 1.0 + paramskappa * (1.0 - exp(-paramsmu * s^2/paramskappa)):
f  := x -> f0(X2S * x):
