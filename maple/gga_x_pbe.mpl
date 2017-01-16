(* type: work_gga_x *)
(* prefix:
  gga_x_pbe_params *params;
 
  assert(p->params != NULL);
  params = (gga_x_pbe_params * )(p->params);
*)

kappa := params_a_kappa:
mu    := params_a_mu:

$include "pbe_x.mpl"

f  := x -> pbe_x(x):
