(* type: work_gga_c *)
(* prefix:
  gga_xc_wb97_params *params;

  assert(p->params != NULL);
  params = (gga_xc_wb97_params * )(p->params);
*)

$define lda_c_pw_params
$include "lda_c_pw.mpl"

$include "lda_x_erf.mpl"

$include "b97.mpl"

f := (rs, z, xt, xs0, xs1) ->
  + f_b97(f_lda_x_erf, 0.004, params_a_c_x, 0, [0, 0, 0, 0, 0],
        rs, z, xs0, xs1)
  + f_b97(f_pw, 0.2, params_a_c_ss, 0.006, params_a_c_ab,
        rs, z, xs0, xs1):

