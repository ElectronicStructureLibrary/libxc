(* type: work_lda *)

params_a_alpha := 1.0:
$include "lda_x.mpl"

beta := rs -> (9.0*Pi/4.0)^(1.0/3.0)/(rs*M_C):
phi  := rs -> 1.0 - 1.5*(sqrt(1.0 + beta(rs)^2)/beta(rs) - arcsinh(beta(rs))/beta(rs)^2)^2:

f    := (rs, z) -> f_lda_x(rs, z)*phi(rs):