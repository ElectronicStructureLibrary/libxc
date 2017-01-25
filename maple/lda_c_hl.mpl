(* type: work_lda *)
(* prefix:
  lda_c_hl_params *params;

  assert(p->params != NULL);
  params = (lda_c_hl_params * )(p->params);
*)

xx   := (k, rs) -> rs/params_a_r[k]:
hl_f := (k, rs) -> -params_a_c[k]*
  ((1.0 + xx(k, rs)^3)*log(1.0 + 1.0/xx(k, rs)) - xx(k, rs)^2 + 0.5*xx(k, rs) - 1.0/3.0):

f := (rs, zeta) -> hl_f(0, rs) + f_zeta(zeta)*(hl_f(1, rs) - hl_f(0, rs)):