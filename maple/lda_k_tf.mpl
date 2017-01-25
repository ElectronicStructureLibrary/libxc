(* type: work_lda *)
(* prefix:
  lda_k_tf_params *params;

  assert(p->params != NULL);
  params = (lda_k_tf_params * )(p->params);
*)

f_zeta_k := z -> 0.5*((1.0 + z)^(5.0/3.0) + (1.0 - z)^(5.0/3.0)):

f := (rs, zeta) -> params_a_ax*f_zeta_k(zeta)/rs^2: