(* type: work_gga_x *)
(* prefix:
  gga_x_pw86_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pw86_params * )(p->params);
*)

$ifdef gga_x_rpw86_params
params_a_aa := 15.0*0.1234:
params_a_bb := 17.33:
params_a_cc := 0.163:
$endif

f0_pw86 := s -> (1.0 + params_a_aa*s^2 + params_a_bb*s^4 + params_a_cc*s^6)^(1.0/15.0):
f_pw86  := x -> f0_pw86(X2S*x):

f      := x -> f_pw86(x):