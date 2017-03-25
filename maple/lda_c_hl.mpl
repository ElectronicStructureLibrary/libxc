(* type: work_lda *)
(* prefix:
  lda_c_hl_params *params;

  assert(p->params != NULL);
  params = (lda_c_hl_params * )(p->params);
*)

$ifdef lda_c_vbh_params
params_a_r := [30, 75.0]:
params_a_c := [0.0252, 0.0127]:
$endif

xx := (k, rs) -> rs/params_a_r[k]:
hl := (k, rs) -> -params_a_c[k]*
  ((1 + xx(k, rs)^3)*log(1 + 1/xx(k, rs)) - xx(k, rs)^2 + 1/2*xx(k, rs) - 1/3):

f_hl := (rs, zeta) -> hl(1, rs) + f_zeta(zeta)*(hl(2, rs) - hl(1, rs)):
f    := (rs, zeta) -> f_hl(rs, zeta):