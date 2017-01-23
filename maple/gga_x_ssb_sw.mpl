(* type: work_gga_x *)
(* prefix:
  gga_x_ssb_sw_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_ssb_sw_params * )(p->params);
*)

f0 := s -> params_a_A 
   + params_a_B*s^2/(1.0 + params_a_C*s^2)
   - params_a_D*s^2/(1.0 + params_a_E*s^4):

f  := x -> f0(X2S*x):