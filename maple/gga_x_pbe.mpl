(* type: work_gga_x *)
(* prefix:
  gga_x_pbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pbe_params * )(p->params);
*)

(* standard PBE *)
$ifdef gga_x_pbe_params
params_a_kappa := 0.8040:
params_a_mu    := 0.2195149727645171:
$endif

(* PBE_SOL *)
$ifdef gga_x_pbe_sol_params
params_a_kappa := 0.8040:
params_a_mu    := MU_GE:
$endif

$ifdef gga_x_pbe_tca_params
params_a_kappa := 1.227:
params_a_mu    := 0.2195149727645171:
$endif

f0_pbe := s -> 1 + params_a_kappa*(1 - params_a_kappa/(params_a_kappa + params_a_mu*s^2)):
f_pbe  := x -> f0_pbe(X2S*x):

f  := x -> f_pbe(x):