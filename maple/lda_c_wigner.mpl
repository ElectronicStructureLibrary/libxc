(* type: work_lda *)
(* prefix:
  lda_c_wigner_params *params;

  assert(p->params != NULL);
  params = (lda_c_wigner_params * )(p->params);
*)

f := (rs, zeta) -> 1*params_a_a/(params_a_b + rs):
