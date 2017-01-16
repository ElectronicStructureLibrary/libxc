(* type: work_mgga_x *)
(* prefix:
  mgga_x_m06l_params *params;
 
  assert(pt->params != NULL);
  params = (mgga_x_m06l_params * )(pt->params);
*)

(* standard PBE *)
kappa := 0.8040:
mu    := 0.2195149727645171:
$include "pbe_x.mpl"

coeff_a := params_a_a:
$include "fw.mpl"

alpha   := 0.00186726:
coeff_d := params_a_d:
$include "gvt4.mpl"

(* there is a factor if 2 in the definition of z, as in Theor. Chem. Account 120, 215 (2008) *)
(* A MINUS was missing in Eq. (7) of the paper *)

f := (rs, x, t, u) -> pbe_x(x)*fw(t, 12) + gtv4(x, 2.0*(t - K_FACTOR_C)):
