(* type: work_gga_x *)
(* prefix:
  gga_x_rpbe_params *params;

  assert(p->params != NULL);
  params = (gga_x_rpbe_params * )(p->params);
*)

$ifdef gga_x_rpbe_params
params_a_rpbe_kappa := KAPPA_PBE:
params_a_rpbe_mu    := MU_PBE:
$endif

f0_rpbe := s -> 1.0 + params_a_rpbe_kappa * (
  1.0 - exp(-params_a_rpbe_mu*s^2/params_a_rpbe_kappa)
):
f_rpbe  := x -> f0_rpbe(X2S*x):

f       := x -> f_rpbe(x):
