(* type: work_mgga_x *)

params_a_mu := MU_GE:
params_a_b  := 0.4:
params_a_c  := 2.14951:
params_a_e  := 1.987:

mgamma      := 0.000023:
ff          := 3:

params_a_kappa := 1/(mgamma/params_a_mu^2 + mgamma/params_a_mu + 1):

$include "tpss_x.mpl"

(* Equation (8) *)

f   := (rs, x, t, u) -> 1 + fx(x, t)*exp(-mgamma*fx(x, t)/params_a_mu)/(1 + fx(x, t)) \
    + (1 - exp(-mgamma*fx(x, t)^2/params_a_mu^2))*(params_a_mu/fx(x, t) - 1):