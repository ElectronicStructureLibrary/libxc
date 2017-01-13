(* type: work_lda *)
(* prefix:
  lda_xc_1d_ehwlrg_params *params;
 
  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);
*)

n := rs -> 1.0/(2.0*rs):

f := (rs, zeta) -> \
 (params_a_a1 + params_a_a2*n(rs) + params_a_a3*n(rs)^2) * n(rs)^params_a_alpha:
