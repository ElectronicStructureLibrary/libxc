(* type: work_lda *)

params_a_alpha := 1:
$include "lda_x.mpl"

beta := rs -> (9*Pi/4)^(1/3)/(rs*M_C):
phi  := rs -> 1 - 1.5*(sqrt(1 + beta(rs)^2)/beta(rs) - arcsinh(beta(rs))/beta(rs)^2)^2:

f    := (rs, z) -> f_lda_x(rs, z)*phi(rs):