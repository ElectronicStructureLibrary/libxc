(* type: work_gga_x *)
(* prefix:
  gga_k_dk_params *params;
 
  assert(p->params != NULL);
  params = (gga_k_dk_params * )(p->params);
*)

f := x -> 
  add(1*params_a_aa[i]*x^(2*(i-1)), i=1..5) /
  add(1*params_a_bb[i]*x^(2*(i-1)), i=1..5):