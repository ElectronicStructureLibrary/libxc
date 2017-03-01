(* type: work_gga_c *)
(* prefix:
  gga_c_bmk_params *params;

  assert(p->params != NULL);
  params = (gga_c_bmk_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$include "b97.mpl"

f := (rs, z, xt, xs_0_, xs_1_) ->
  f_b97(f_pw, 0.2, params_a_c_ss, 0.006, params_a_c_ab,
        rs, z, xs_0_, xs_1_):
