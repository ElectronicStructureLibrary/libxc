(* type: work_lda *)
(* prefix:
  lda_c_2d_prm_params *params;

  assert(p->params != NULL);
  params = (lda_c_2d_prm_params * )(p->params);

  assert(params->N > 1.0);
*)

prm_q := 3.9274: (* 2.258 *)

(* Equation (4) *)
beta := rs -> prm_q/(sqrt(Pi)*rs):

phi  := rs -> beta(rs)/(beta(rs) + sqrt(Pi)/2.0):

f0 := rs ->
  + sqrt(Pi)*beta(rs)*(phi(rs) - 1.0)/(2.0*sqrt(2.0 + params_a_c)) (* original version has (phi-1)^2 *)
  + phi(rs)*(phi(rs) - 1.0)/(2.0 + params_a_c)
  + sqrt(Pi)*phi(rs)*phi(rs)/(4.0*beta(rs)*(2.0 + params_a_c)^1.5)
  + sqrt(Pi)*beta(rs)*(phi(rs) - 1.0)/sqrt(1.0 + params_a_c)
  + phi(rs)/(1.0 + params_a_c):

f  := (rs, z) -> f0(rs)*Pi/(2.0*prm_q*prm_q):
