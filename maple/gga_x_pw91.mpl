(* type: work_gga_x *)
(* prefix:
  gga_x_pw91_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw91_params * )(p->params);
*)

pw91_num := s -> (params_a_c + params_a_d*exp(-params_a_alpha*s^2))*s^2
         - params_a_f*s^params_a_expo:
pw91_den := s -> 1.0 + s*params_a_a*arcsinh(params_a_b*s) + params_a_f*s^params_a_expo:

f  := x -> 1.0 + pw91_num(X2S*x)/pw91_den(X2S*x):
