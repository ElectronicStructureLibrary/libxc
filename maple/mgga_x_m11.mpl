(* type: work_mgga_c *)
(* prefix:
  mgga_x_m11_params *params;

  assert(p->params != NULL);
  params = (mgga_x_m11_params * ) (p->params);
*)

$include "mgga_x_m08.mpl"
$include "lda_x_erf.mpl"

f_spin := (rs, z, x, t) ->
  + lda_x_ax*(1.0 + z)^(4.0/3.0)/rs
  * attenuation_erf(a_cnst*rs/(1.0 + z)^(1/3))
  * m08_f(params_a_a, params_a_b, x, t):

f_m11 := (rs, z, xt, xs0, xs1, ts0, ts1) ->
  f_spin(rs, z, xs0, ts0) + f_spin(rs, -z, xs1, ts1):

f := (rs, z, xt, xs0, xs1, ts0, ts1, us0, us1) ->
  f_m11(rs, z, xt, xs0, xs1, ts0, ts1):
