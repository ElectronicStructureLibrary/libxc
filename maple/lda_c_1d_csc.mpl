(* type: work_lda *)
(* prefix:
  lda_c_1d_csc_params *params;

  assert(p->params != NULL);
  params = (lda_c_1d_csc_params * )(p->params);
*)

(* factor of 2 is the conversion from Ry to Hartree *)
f_aux := (a, rs) -> -(rs + a[5]*rs^2)*log(1+ a[8]*rs + a[9]*rs^a[10])
  / (2*(a[1] + a[2]*rs + a[3]*rs^a[6] + a[4]*rs^a[7])):

f := (rs, z) ->
  f_aux(params_a_para, rs) + (f_aux(params_a_ferro, rs) - f_aux(params_a_para, rs))*z^2:
