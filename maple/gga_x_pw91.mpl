(* type: work_gga_x *)
(* prefix:
  gga_x_pw91_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw91_params * )(p->params);
*)

$ifdef gga_x_pw91_params
params_a_a     :=   0.19645:
params_a_b     :=   7.7956:
params_a_c     :=   0.2743:
params_a_d     :=  -0.1508:
params_a_f     :=   0.004:
params_a_alpha := 100:
params_a_expo  :=   4:
$endif

pw91_num := s -> (params_a_c + params_a_d*exp(-params_a_alpha*s^2))*s^2
         - params_a_f*s^params_a_expo:
pw91_den := s -> 1 + s*params_a_a*arcsinh(params_a_b*s) + params_a_f*s^params_a_expo:

f_pw91  := x -> 1 + pw91_num(X2S*x)/pw91_den(X2S*x):
f       := x -> f_pw91(x):