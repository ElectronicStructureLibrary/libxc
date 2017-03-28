(* type: work_lda *)
(* prefix:
  lda_c_wigner_params *params;

  assert(p->params != NULL);
  params = (lda_c_wigner_params * )(p->params);
*)

f := (rs, z) -> (1 - z^2)*params_a_a/(params_a_b + rs):
