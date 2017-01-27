(* type: work_lda *)
(* prefix:
  lda_x_params *params;

  assert(p->params != NULL);
  params = (lda_x_params * )(p->params);
*)

ax := -params_a_alpha*RS_FACTOR*X_FACTOR_C/2.0^(4.0/3.0):

f_lda_x := (rs, z) -> ax*((1.0 + z)^(4.0/3.0) + (1.0 - z)^(4.0/3.0))/rs:
f       := (rs, z) -> f_lda_x(rs, z):