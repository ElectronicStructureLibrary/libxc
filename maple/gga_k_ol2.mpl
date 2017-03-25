(* type: work_gga_x *)
(* prefix:
  gga_k_ol2_params *params;
 
  assert(p->params != NULL);
  params = (gga_k_ol2_params * )(p->params);
*)

f := x -> 
  + params_a_aa
  + params_a_bb*x^2/72.0
  + params_a_cc*x/(2^(1/3) + 4*x):